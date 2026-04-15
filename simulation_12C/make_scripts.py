import sys
import os

sys.path.append("../common")
from atomx_simulation_script_writter import atomx_simulation_script_writter

path_to_data = "data"
path_to_figures = "figures"
path_to_input = "input"

if __name__ == "__main__":
    #print(atomx_simulation_script_writter.__doc__)

    a = atomx_simulation_script_writter("ATOMX_12C")
    a.set_path(path_to_data, path_to_input, path_to_figures)
    a.make_project_configuration()
    a.make_geant4_vis_file()

    threshold = 7.65
    detector = a.make_detector_file(reaction_z=[150,160], step_limit=10)
    reaction = a.make_eventgen_file_beam(beam_particle="1H", energy=12)
    a.add_eventgen_elastic(reaction, beam_particle="1H", target="12C", heavy="12C", light="1H", ex_energy_heavy=threshold)
    a.add_eventgen_decay  (reaction, threshold=threshold, daughters=["4He","4He","4He"], ex_energies=[0,0,0], branching_ratio=1, life_time=0, shoot=[1,1,1])

    a.make_script(detector.name, reaction.name)
