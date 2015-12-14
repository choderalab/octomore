#!/usr/bin/env python

import openpathsampling as paths
import numpy as np
import mdtraj as md
import math

Abl = {
    'name' : "Abl",
    'file' : "sim-snippets/dozen_frames_abl.xtc",
    'pdb' : "sim-snippets/abl_ref.pdb",
    'DFG' : [2257,2255,2265,2270]
}

class Unimore(object):
    """If Octomore is about finding 8+ possible CVs, this is for 1 of them"""
    def build_engine(self):
        """Builds the OpenMMEngine."""

        from simtk.openmm import app
        import simtk.openmm as mm
        from simtk import unit
        kinase = self.kinase

        #### SET UP THE OPENMM STUFF #######################################
        # taken from builder.openmm.org
        forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

        system = forcefield.createSystem(app.PDBFile(kinase['pdb']).topology,
                                         nonbondedMethod=app.NoCutoff, 
                                         constraints=app.HBonds)

        integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                           1.0/unit.picoseconds, 
                                           2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)


        #### OPENPATHSAMPLING-SPECIFIC SETUP ###############################
        options = {'nsteps_per_frame' : 50, 'n_frames_max' : 1500}
        template = paths.tools.snapshot_from_pdb(kinase['pdb'])

        engine = paths.OpenMMEngine(
                template=template,
                system=system,
                integrator=integrator,
                options=options
        )

        return engine

    def get_initial_frame(self, frame_num, file_name, pdb):
        """Pulls an initial frame out of an xtc file"""
        traj = paths.trajectory_from_mdtraj(md.load(file_name, top=pdb))
        return traj[frame_num]

    def __init__(self, output_file=None, kinase=Abl):
        self.cv = None # override in subclasses
        self.kinase = kinase
        self.engine = self.build_engine()
        self.engine.current_snapshot = self.get_initial_frame(
            frame_num=1,
            file_name=kinase['file'],
            pdb=kinase['pdb']
        )
        # REALLY? This is ridiculous, and illustrates the problem with the
        # topology approach
        self.engine.current_snapshot.configuration.topology = self.engine.topology

        self.dfg = paths.CV_MDTraj_Function(name="DFG", 
                                            f=md.compute_dihedrals,
                                            indices=[kinase['DFG']])
        self.DFG_out = paths.CVRangeVolumePeriodic(
            self.dfg,
            lambda_min=-4.0, lambda_max=-2.6,
            period_min=-math.pi, period_max=math.pi
        )
        self.DFG_in = paths.CVRangeVolumePeriodic(
            self.dfg,
            lambda_min=0.0, lambda_max=1.0,
            period_min=-math.pi, period_max=math.pi
        )
        if output_file is not None:
            self.storage = paths.storage.Storage(
                filename=output_file,
                mode="w",
                template=self.engine.template
            )
        else:
            self.storage = None

    def ratcheter(self, interfaces, direction="out_in"):
        if direction == "in_out":
            stateA = self.DFG_in
            stateB = self.DFG_out
        elif direction == "out_in":
            stateA = self.DFG_out
            stateB = self.DFG_in
        else:
            raise RuntimeError("direction is " + str(direction) + 
                               ": must be 'in_out' or 'out_in'.")
        transition = paths.TISTransition(stateA, stateB, interfaces, self.cv)
        ratcheter = paths.FullBootstrapping(
            transition=transition,
            snapshot=self.engine.current_snapshot,
            storage=self.storage,
            engine=self.engine
        )
        return ratcheter


class DFG(Unimore):
    """As an example CV, we implement the DFG (which is already prepared)"""
    def __init__(self, output_file=None, kinase=Abl):
        super(DFG, self).__init__(output_file, kinase)
        self.cv = self.dfg

