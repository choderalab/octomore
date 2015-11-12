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

        # taken from builder.openmm.org (working in vacuum here)
        forcefield = app.ForceField('amber99sbildn.xml')

        system = forcefield.createSystem(app.PDBFile(Abl_topol).topology,
                                         nonbondedMethod=app.PME, 
                                         nonbondedCutoff=1.0*unit.nanometers,
                                         constraints=app.HBonds,
                                         rigidWater=True, 
                                         ewaldErrorTolerance=0.0005)
        integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                           1.0/unit.picoseconds, 
                                           2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)

        options = { 'nsteps_per_frame' : 10}

        engine = paths.OpenMMEngine(
                template=template,
                system=system,
                integrator=integrator,
                options=options
        )
        return engine

    def get_initial_frame(self, frame_num, file_name, pdb):
        """Pulls an initial frame out of an xtc file"""
        traj = paths.Trajectory.from_mdtraj(md.load(file_name, top=pdb))
        return traj[frame_num]

    def __init__(self, output=None, kinase=Abl):
        self.cv = None # override in subclasses
        self.engine = self.build_engine()
        self.engine.current_snapshot = get_initial_frame(
            frame_num=0,
            file_name=kinase['file'],
            top=kinase['pdb']
        )
        self.dfg = paths.CV_MDTraj_Function(name="DFG", 
                                            f=md.compute_dihedrals,
                                            indices=kinase['DFG'])
        self.DFG_in = paths.CVRangeVolumePeriodic(
            dfg,
            lambda_min=-3.1, lambda_max=0.0,
            period_min=-math.pi, period_max=math.pi
        )
        self.DFG_out = paths.CVRangeVolumePeriodic(
            dfg
            lambda_min=2.1, lambda_max=3.0,
            period_min=-math.pi, period_max=math.pi
        )
        self.storage = paths.storage.Storage(
            filename=output,
            mode="w",
            template=self.engine.current_snapshot,
        )

    def ratcheter(self, interfaces, direction="in_out"):
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
        ratcheter = paths.FullBoostrapping(
            transition=transition,
            snapshot=self.engine.current_snapshot,
            storage=self.storage,
            engine=self.engine
        )
        return ratcheter


class DFG(Unimore):
    def __init__(self, kinase=Abl):
        super(DFG, self).__init__(kinase)
        self.cv = self.dfg

