import os
import math
import numpy as np

class atomx_simulation_script_writter:
    """
# List of methods
- def set_path(self, input_path="input", data_path="data", figures_path="figures", exist_ok=True)
- def make_script(self, detector_file_name, reaction_file_name, macro_for_sim="#", macro_for_ana="#")
- def make_detector_file(self, detector_file_name="", make_gas=True, make_si_side=True, make_si_downstream=True, make_si_corner=True, make_si_bottom=True, make_csi=True):
- def make_eventgen_file_beam(self, beam_particle, energy, eventgen_name="", z_emission=-1000, sigma_energy=0, sigma_theta=0, sigma_x=0, sigma_y=0):
- def make_eventgen_file_isotropic(self, particle_name, energy1, energy2, angle1, angle2, z0=0, sz=0, eventgen_name)
- def make_cs_file_gaus(self, mean, sigma)
- def make_cs_file_from_data(self, cs_file_name, angle1, angle2)
- def add_eventgen_elastic(self, reaction_file, beam_particle, target, heavy, light, ex_energy_light, ex_energy_heavy, cs_file_name = "flat.txt"):
- def add_eventgen_decay  (self, reaction_file, threshold=7.65, daughters=["4He","4He","4He"], ex_energies=[0,0,0], branching_ratio=1, life_time=0, shoot=[1,1,1]):
- def make_geant4_vis_file(self)
- def make_project_configuration(self)
"""

    def __init__(self, name):
        self.input_path = "input"
        self.data_path = "data"   
        self.figures_path = "figures"
        self.project_name = name

    def set_path(self, data_path="data", input_path="input", figures_path="figures"):
        self.input_path = input_path  
        self.data_path = data_path   
        self.figures_path = figures_path

        for key, link_path in [['data',data_path],['input',input_path],['figures',figures_path]]:
            key_exist = os.path.exists(key)
            key_isdir = os.path.isdir(key)
            key_islink = os.path.islink(key)
            key_wblink = (link_path.find('/')>=0)
            if key_exist and key_isdir and key_islink==False:
                print(f'!! "{key}" is a directory! Remove it manually and run the script again (if you have to)')
            else:
                if key_exist and key_islink:
                    print(f'-- removing link {key}')
                    os.remove(key)
                if key_wblink:
                    print(f'-- creating symbolic link {link_path} to {key}')
                    os.symlink(link_path, key)
                else:
                    print(f'-- creating directory {data_path}')
                    os.makedirs(data_path, exist_ok=True)

    def make_detector_file(self, detector_file_name="", reaction_z=[], step_limit=0.5, make_gas=True, make_si_side=True, make_si_downstream=True, make_si_corner=True, make_si_bottom=True, make_csi=True):
        box_size = [480, 450, 480]
        gas_size = [450, 400, 450]
        pad_size = [400, 2,   400]
        mms_size = [400, 20,  400]
        gas_name = ["He", "CO2"]
        gas_perc = [97, 3]
        if len(reaction_z)==0: reaction_z = [0, 0]
        #pcb_frame_size = [258, 300, 290]
        larger_frame_size = [280, 300, 320]
        x6_size = [45.20,93.10,2.40]
        side_plane_z_center = 0
        downst_plane_z = larger_frame_size[2]/2
        bottom_side_x = x6_size[0] + 0.5*x6_size[1]
        bottom_z_end = 130
        bottom_y = -larger_frame_size[1]*0.45
        pos_xyz_corner_list = [
            [+110, 0, +120, -135],
            [-110, 0, +120,  135],
            [-110, 0, -120,  45 ],
            [+110, 0, -120, -45 ],
        ]

        if make_si_side:
            det_par_list = []
            det_type = "X6"
            pos_xyz_side = [larger_frame_size[0]/2, 0, side_plane_z_center]
            for side in range(2):
                for layer in range(2):
                    for row in range(4):
                        pos_xyz = [pos_xyz_side[0],pos_xyz_side[1],pos_xyz_side[2]]
                        pos_xyz[2] = pos_xyz[2] + (row-1.5)*x6_size[0]
                        pos_xyz[1] = pos_xyz[1] + 0.5*x6_size[1]
                        if layer==1: pos_xyz[1] = pos_xyz[1] - 1.0*x6_size[1]
                        rot_xyz = [180,90,0]
                        if layer==0: rot_xyz[1] = rot_xyz[1] + 180
                        if layer==1: rot_xyz[0] = 0
                        if side==1:
                            pos_xyz[0] = -pos_xyz[0]
                            rot_xyz[1] = rot_xyz[1] + 180
                        name = f"side_{side}_{layer}_{row}"
                        csi_side = 1 if make_csi else 0
                        det_par_list.append([name,det_type,pos_xyz,rot_xyz,csi_side])

        if make_si_downstream:
            pos_xyz_downst = [0, 0, downst_plane_z]
            for layer in range(2):
                for row in range(3):
                    pos_xyz = [pos_xyz_downst[0],pos_xyz_downst[1],pos_xyz_downst[2]]
                    pos_xyz[0] = pos_xyz[0] + (row-1)*x6_size[0]
                    pos_xyz[1] = pos_xyz[1] + 0.5*x6_size[1]
                    if layer==1: pos_xyz[1] = pos_xyz[1] - 1.0*x6_size[1]
                    rot_xyz = [180,0,0]
                    if layer==0: rot_xyz[1] = rot_xyz[1] + 180
                    if layer==1: rot_xyz[0] = 0
                    name = f"downstream_{layer}_{row}"
                    csi_side = 1 if make_csi else 0
                    det_par_list.append([name,det_type,pos_xyz,rot_xyz,csi_side])

        if make_si_corner:
            for i in range(len(pos_xyz_corner_list)):
                for layer in range(2):
                    pos_xyz = [pos_xyz_corner_list[i][0],pos_xyz_corner_list[i][1],pos_xyz_corner_list[i][2]]
                    rot_xyz = [0,pos_xyz_corner_list[i][3],0]
                    pos_xyz[1] = pos_xyz[1] + 0.5*x6_size[1]
                    if layer==1: pos_xyz[1] = pos_xyz[1] - 1.0*x6_size[1]
                    if layer==0: rot_xyz[0] = 180
                    if layer==1:
                        rot_xyz[0] = 0
                        rot_xyz[1] = rot_xyz[1] + 180
                    name = f"corner_{i}_{layer}"
                    csi_side = 1 if make_csi else 0
                    det_par_list.append([name,det_type,pos_xyz,rot_xyz,csi_side])

        if make_si_bottom:
            pos_xyz = [-0.5*x6_size[0], bottom_y, bottom_z_end-x6_size[1]*0.42]
            rot_xyz = [90,0,0]
            name = f"bottomds_0"
            csi_side = 1 if make_csi else 0
            det_par_list.append([name,det_type,pos_xyz,rot_xyz,csi_side])
            pos_xyz = [+0.5*x6_size[0], bottom_y, bottom_z_end-x6_size[1]*0.42]
            rot_xyz = [90,0,0]
            name = f"bottomds_1"
            csi_side = 1 if make_csi else 0
            det_par_list.append([name,det_type,pos_xyz,rot_xyz,csi_side])
            pos_xyz_bottom_side = [0, bottom_y, bottom_z_end]
            for side in range(2):
                for row in range(7):
                    pos_xyz = [pos_xyz_bottom_side[0],pos_xyz_bottom_side[1],pos_xyz_bottom_side[2]]
                    pos_xyz[2] = pos_xyz[2] - (row+0.5)*x6_size[0]
                    if side==0: pos_xyz[0] = pos_xyz[0] - bottom_side_x
                    if side==1: pos_xyz[0] = pos_xyz[0] + bottom_side_x
                    rot_xyz = [0,0,0]
                    if side==0: rot_xyz = [90,90 ,0]
                    if side==1: rot_xyz = [90,270,0]
                    name = f"bottomside_{side}_{row}"
                    csi_side = 1 if make_csi else 0
                    det_par_list.append([name,det_type,pos_xyz,rot_xyz,csi_side])

        if len(detector_file_name)==0:
            detector_file_name = f"{self.input_path}/{self.project_name}.detector"
        print("-- creating", detector_file_name)
        detector_file = open(detector_file_name,"w")
        print(f"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detector file {detector_file_name} ({self.project_name})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%""", file=detector_file)
        if make_gas:
            print(f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gas
ATOMX
    TPC_Pos              = 0 0 0 mm
    TPC_Chamber_Size     = {box_size[0]} {box_size[1]} {box_size[2]} mm
    TPC_Gas_Size         = {gas_size[0]} {gas_size[1]} {gas_size[2]} mm
    TPC_Pad_Size         = {pad_size[0]} {pad_size[1]} {pad_size[2]} mm
    TPC_MMS_Size         = {mms_size[0]} {mms_size[1]} {mms_size[2]} mm
    TPC_Mylar_Rmax       = 3.5 cm
    TPC_Mylar_Thickness  = 7 micrometer
    Gas_EProductionCut   = 1000 mm
    Gas_StepLimit        = {step_limit} mm
    Gas_Material         = {gas_name[0]} {gas_name[1]}
    Gas_Fraction         = {gas_perc[0]} {gas_perc[1]} 
    Gas_Temperature      = 295 kelvin
    Gas_Pressure         = 0.1 bar
    Gas_ReactionZ        = {reaction_z[0]} {reaction_z[1]} mm
    Gas_ReactionBox_Size = 50 mm""", file=detector_file)

        for det_par in det_par_list:
            name = det_par[0]
            det_type = det_par[1]
            pos_xyz = det_par[2]
            rot_xyz = det_par[3]
            csi_side = det_par[4]
            mother_volume_parameter = "MotherVolume = ATOMX_LV_Gas" if make_gas else "%MotherVolume"
            print(f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% X6 ({name})
STARK
    {mother_volume_parameter}
    Type = {det_type}
    POS = {pos_xyz[0]} {pos_xyz[1]} {pos_xyz[2]} mm
    RotateXYZ = {rot_xyz[0]} {rot_xyz[1]} {rot_xyz[2]} deg
    CsI = {csi_side}""", file=detector_file)

        return detector_file

    def make_cs_file_gaus(self, mean, sigma):
        cs_file_name = f"{self.input_path}/gaus_{mean}_{sigma}.txt"
        print("-- creating", cs_file_name)
        cs_file = open(cs_file_name,'w')
        for i in range(0, 1800):
            x = i/10
            value = np.exp(-(x - mean)**2 / (2 * sigma**2))
            cs_file.write(f"{x}\t{value}\n")
        cs_file.write(f"180 0")
        return cs_file_name

    def make_cs_file_from_data(self, cs_file_name, angle1, angle2):
        cs_file = open(cs_file_name)
        file0 = open(cs_file_name)
        cs_file_name1 = f'{cs_file_name[:cs_file_name.rfind(".")]}_{angle1}_{angle2}{cs_file_name[cs_file_name.rfind("."):]}'
        print("-- creating", cs_file_name1)
        file1 = open(cs_file_name1,'w')
        for line in file0:
            theta, value = line.split()
            if float(theta)>=angle1 and float(theta)<angle2:
                print(theta,value,file=file1)
            else:
                print(theta,0,file=file1)
        return file1

    def make_eventgen_file_beam(self, beam_particle, energy, eventgen_name="", z_emission=-1000, sigma_energy=0, sigma_theta=0, sigma_x=0, sigma_y=0):
        if len(eventgen_name)==0: eventgen_name = f"{self.project_name}"
        if eventgen_name.endswith(".reaction")==False: eventgen_name = f"{eventgen_name}.reaction"
        if eventgen_name.startswith(f"{self.input_path}/")==False: eventgen_name = f"{self.input_path}/{eventgen_name}"
        reaction_file = open(eventgen_name,'w')
        print("-- creating", eventgen_name)
        content = f"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Event generator for {eventgen_name} ({self.project_name})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beam
    Particle= {beam_particle}
    Energy= {energy} MeV
    SigmaEnergy= {sigma_energy} MeV
    SigmaThetaX= {sigma_theta} deg
    SigmaPhiY= 0 deg
    SigmaX= {sigma_x} mm
    SigmaY= {sigma_y} mm
    MeanThetaX= 0. deg
    MeanPhiY= 0 deg
    MeanX= 0 mm
    MeanY= 0 mm
    ZEmission = {z_emission} mm"""
        print(content, file=reaction_file)
        return reaction_file

    def add_eventgen_elastic(self, reaction_file, beam_particle, target, heavy, light, ex_energy_light=0, ex_energy_heavy=0, cs_file_name = "flat.txt"):
        content = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TwoBodyReaction
    Beam= {beam_particle}
    Target= {target}
    Light= {light}
    Heavy= {heavy}
    ExcitationEnergyLight= {ex_energy_light} MeV
    ExcitationEnergyHeavy= {ex_energy_heavy} MeV
    CrossSectionPath= {cs_file_name} CS
    ShootLight= 1
    ShootHeavy= 1"""
        print(content, file=reaction_file)
        return reaction_file

    def add_eventgen_decay(self, reaction_file, threshold=7.65, daughters=["4He","4He","4He"], ex_energies=[0,0,0], branching_ratio=1, life_time=0, shoot=[1,1,1]):
        ex_energies = [str(x) for x in ex_energies]
        shoot = [str(x) for x in shoot]
        print(f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Decay 12C
    Threshold = {threshold} MeV
    Daughter = {' '.join(daughters)}
    ExcitationEnergy = {' '.join(ex_energies)} MeV
    BranchingRatio = {branching_ratio}
    LifeTime = {life_time} ns
    Shoot = 1 1 1
""", file=reaction_file)

    def make_eventgen_file_isotropic(self, particle_name, energy1, energy2, angle1, angle2, dist="flat", x0=0, y0=0, z0=0, sx=0, sy=0, sz=0, ex_energy=0., eventgen_name=""):
        if len(eventgen_name)==0: eventgen_name = f"{self.project_name}"
        if eventgen_name.endswith(".reaction")==False: eventgen_name = f"{eventgen_name}.reaction"
        if eventgen_name.startswith(f"{self.input_path}/")==False: eventgen_name = f"{self.input_path}/{eventgen_name}"
        reaction_file = open(eventgen_name,'w')
        print("-- creating", eventgen_name)
        content = f"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isotropic event generator {eventgen_name} ({self.project_name})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Isotropic
    Particle                      = {particle_name}
    EnergyDistribution            = {dist}
    EnergyLow                     = {energy1} MeV
    EnergyHigh                    = {energy2} MeV
    HalfOpenAngleMin              = {angle1} deg
    HalfOpenAngleMax              = {angle2} deg
    x0                            = {x0} mm
    y0                            = {y0} mm
    z0                            = {z0} mm
    ExcitationEnergy              = {ex_energy} MeV
    SigmaX                        = {sx} mm
    SigmaY                        = {sy} mm
    SigmaZ                        = {sz} mm
    Multiplicity                  = 1
    %Direction                     =
    %EnergyDistributionHist        =
    %SourceProfile                 ="""
        print(content, file=reaction_file)
        return reaction_file

    def make_script(self, detector_file_name, reaction_file_name, run_name="", nbeams=1000, macros1="#", macros2="#", macro_for_sim="#", macro_for_ana="#", run_make_class_tree=False):
        if len(run_name)==0: run_name = self.project_name
        if type(macros1)==list: macros1 = '\n'.join(macros1)
        if type(macros2)==list: macros2 = '\n'.join(macros2)
        batch_file_name = f"batch_{run_name}.sh"
        viewer_file_name = f"viewer_{run_name}.sh"
        batch_file = open(batch_file_name,'w')
        viewer_file = open(viewer_file_name,'w')
        os.chmod(batch_file_name, 0o755)
        os.chmod(viewer_file_name, 0o755)
        print("-- creating", batch_file_name)
        print("-- creating", viewer_file_name)
        content = f"#!/bin/bash\n"
        print(content, file=batch_file)
        print(content, file=viewer_file)
        sim_file_name = f"{run_name}.sim.root"
        ana_file_name = f"{run_name}.ana.root"
        vwr_file_name = f"{run_name}.viewer.root"
        content_b = f"""tee > {self.input_path}/geant4_batch.mac <<EOF
/run/beamOn {nbeams}
EOF

npsimulation --record-track -D {detector_file_name} -E {reaction_file_name} -O {sim_file_name} -B {self.input_path}/geant4_batch.mac
npanalysis -T {self.data_path}/{sim_file_name} SimulatedTree -O {ana_file_name}"""
        content_v = f"""npsimulation -D {detector_file_name} -E {reaction_file_name} -O {vwr_file_name} -M {self.input_path}/geant4_vis.mac"""
        print(content_b, file=batch_file)
        print(content_v, file=viewer_file)
        if run_make_class_tree:
            print(f"""
echo root -l -q -b -e '((TTree*)(new TFile("{self.data_path}/{sim_file_name}"))->Get("SimulatedTree"))->MakeClass();'
root -l -q -b -e '((TTree*)(new TFile("{self.data_path}/{sim_file_name}"))->Get("SimulatedTree"))->MakeClass();'""",file=batch_file)
        print(f"""
{macros1}
{macro_for_sim} -- \\"{self.data_path}/{sim_file_name}\\"
{macro_for_ana} -- \\"{self.data_path}/{ana_file_name}\\"
{macros2}""",file=batch_file)
        return batch_file, viewer_file

    def make_geant4_vis_file(self):
        geant4_vis_file_name = f"{self.input_path}/geant4_vis.mac"
        print("-- creating", geant4_vis_file_name)
        geant4_vis_file = open(geant4_vis_file_name,'w')
        print(f"""/vis/open OGL 600x600-0+0
/vis/scene/add/axes 0 0 0 100 mm
/vis/drawVolume
/vis/viewer/flush
/vis/viewer/set/background white
/vis/viewer/set/lightsVector -1 0 0
/vis/viewer/set/style surface
/vis/viewer/set/auxiliaryEdge true
/vis/viewer/set/lineSegmentsPerCircle 100
/vis/scene/add/trajectories smooth
/vis/modeling/trajectories/create/drawByCharge
/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true
/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2
/vis/scene/endOfEventAction accumulate
/vis/geometry/set/visibility logicWorld 0 false
/vis/viewer/set/viewpointThetaPhi 100 160
/vis/viewer/set/autoRefresh true
/tracking/verbose 1
/vis/verbose warnings
/run/setCutForAGivenParticle e- 100. mm
/run/beamOn 1""", file=geant4_vis_file)

    def make_project_configuration(self):
        proj_config_file_name = "project.config"
        print("-- creating", proj_config_file_name)
        proj_config_file = open(proj_config_file_name,'w')
        print(f"""Project {self.project_name}
    AnalysisOutput   = {self.data_path}
    SimulationOutput = {self.data_path}
    %CalibrationOutput
    %EfficiencyOutput
    %EnergyLoss
    %Cuts""",file=proj_config_file)

        physics_list_file_name = "PhysicsListOption.txt"
        print("-- creating", physics_list_file_name)
        physics_list_file = open(physics_list_file_name,'w')
        print(f"""EmPhysicsList                   Option4 % INCLXX_EM Option1 Option2 Option3 Option4 Standard Livermore Penelope

DefaultCutOff                   1 % mm
IonBinaryCascadePhysics         0
NPIonInelasticPhysics           0
EmExtraPhysics                  0
HadronElasticPhysics            0
StoppingPhysics                 0
OpticalPhysics                  0
DriftElectronPhysics            0

HadronPhysicsQGSP_BIC_HP        0
HadronPhysicsQGSP_BERT_HP       0
HadronPhysicsINCLXX             0
HadronPhysicsQGSP_INCLXX_HP     0
HadronPhysicsQGSP_INCLXX        0
HadronPhysicsFTFP_INCLXX_HP     0
HadronPhysicsFTFP_INCLXX        0
Decay                           0

IonGasModels                    0
pai                             0
pai_photon                      0
Menate_R                        0
NeutronHP                       0
Shielding                       0
Shielding_HPT                   0
ShieldingLEND                   0
%LevelData                       #Z A File""", file=physics_list_file)

        cleaner_file_name = "clean.sh"
        print("-- creating", cleaner_file_name)
        cleaner_file = open(cleaner_file_name,'w')
        print(f"""set -x
rm -f batch_*.sh
rm -f viewer_*.sh
rm -f PhysicsListOption.txt
rm -f project.config
rm -f {self.input_path}/*.reaction
rm -f {self.input_path}/*.detector
rm -f {self.input_path}/geant4*.mac
rm -f clean.sh""", file=cleaner_file)
