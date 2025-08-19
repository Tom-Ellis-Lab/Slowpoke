
# Golden Gate Cloning protocol
# Written by Fankang Meng and Koray Malci, Imperial College London
# Updated by Alicia Da Silva and Henri Galez for the Flex robot, Institut Pasteur

from opentrons import protocol_api, types
import time
import math

metadata = {
    'protocolName': 'GG Cloning',
    'description': 'GG & transformation & plating using a Flex robot for genetic toolkits.'}

requirements = {"robotType": "Flex", "apiLevel": "2.21"}

num_rxns = len(combinations_to_make)
volume_buffer = 1
volume_enzyme = 1
volume_reaction = 10
nb_parts = 4
volume_inputDNA = 1
volume_waterbuffer_per_reaction = volume_reaction - volume_buffer - nb_parts * volume_inputDNA
volume_master_mix = volume_enzyme + volume_waterbuffer_per_reaction
volume_competent_cells = 50

temp_reaction = 4
temp_reagent = 4

def run(protocol: protocol_api.ProtocolContext):
    # Trash needs to be specified with Flex
    trash = protocol.load_trash_bin("A3")

    # Load in 1 10ul tiprack and 2 300ul tipracks
    tr_300 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B3')
    tr_20 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'C3')

    # Load in pipettes
    p50_single = protocol.load_instrument('flex_1channel_50', 'right', tip_racks=[tr_20,tr_300])

    # ========== this block should be activated if a thermocycler module is used in the protocol ==========
    # Load in Bio-Rad 96 Well Plate on Thermocycler Module for GG Assembly, transformation, and outgrowth.
    #tc_mod = protocol.load_module('Thermocycler Module')
    #reaction_plate = tc_mod.load_labware('biorad_96_wellplate_200ul_pcr')
    #tc_mod.open_lid()
    #tc_mod.set_block_temperature(4)
    
    # ========== this block should be activated if an additional temperature module is used instead of a thermocycler module  ==========
    temp_mod_reaction = protocol.load_module('temperature module gen2', 'A1')
    temp_mod_reaction.set_temperature(celsius=temp_reaction)
    reaction_plate = temp_mod_reaction.load_labware('biorad_96_wellplate_200ul_pcr', 'A1')



    # Load in Water & Enzymer+ Buffer, wash water, dilution water, and wash trough (USA Scientific 12 Well Reservoir 22ml)
    temp_mod = protocol.load_module('temperature module gen2', 'D3')
    temp_mod.set_temperature(celsius=temp_reagent)
    trough = temp_mod.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap', 'D3')
    buffer_H2O_mix = trough.wells()[0]  # Well A1
    well_enzyme = trough.wells()[1]  # Well B1
    dilution_water = trough.wells()[2]  # Well C1
    competent_cell = trough.wells()[3]  # Well D1
    liquid_waste = trough.wells()[4]  # Well A2

    # Load in Input DNA Plate
    dna_plate_dict = {}
    plate_name = list(dna_plate_map_dict.keys())
    dna_plate_dict[plate_name[0]] = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'C2', 'Input DNA Plate')
    dna_plate_dict[plate_name[1]] = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 'D2', 'Input DNA Plate2')

    # Load in Agar plate
    agar_plate = protocol.load_labware('corning_6_wellplate_16.8ml_flat', 'C1', 'Agar Plate')


    # This function checks the existence of DNA parts and returns for well location of the parts
    def find_dna(name, dna_plate_map_dict, dna_plate_dict):
        """Return a well containing the named DNA."""
        well_num=None
        rows = ['A','B','C','D','E','F','G','H','I','J']
        columns = [1,2,3,4,5,6,7,8,9,10,11,12]
        for plate_name, plate_map in dna_plate_map_dict.items():
            for i, row in enumerate(plate_map):
                for j, dna_name in enumerate(row):
                    if dna_name == name:
                        plate_part = plate_name
                        if plate_name == 'fixed_toolkit_map':
                            #well_num = 8 * j + i
                            well_num = rows[i]+str(columns[j])
                        elif plate_name == 'custom_parts_map':
                            #well_num = 4 * j + i
                            well_num = rows[i]+str(columns[j])
                         # Check if well_num is assigned a value
        if well_num is not None:
            return dna_plate_dict[plate_part].wells_by_name()[well_num]
            #return [plate_part,well_num]

        # Handle the case where the name is not found
        print(f"DNA with name '{name}' not found.")
        return None  # or raise an exception, depending on your requirements   
                        #return dna_plate_dict[plate_name].wells()[well_num]
        raise ValueError("Could not find dna piece named \"{0}\"".format(name))

    # This function checks if the DNA parts exist in the DNA plates and returns for well location of output DNA combinations
    def find_combination(name, combinations_to_make):
        """Return a well containing the named combination."""
        for i, combination in enumerate(combinations_to_make):
            if combination["name"] == name:
                return reaction_plate.wells()[i]
        raise ValueError("Could not find combination \"{0}\".".format(name))

    combinations_by_part = {}
    for i in combinations_to_make:
        name = i["name"]
        for j in i["parts"]:
            if j in combinations_by_part.keys():
                combinations_by_part[j].append(name)
            else:
                combinations_by_part[j] = [name]

    # This section will take the GG buffer and water into the designation wells
    protocol.pause('Temperature modules ready!')
    p50_single.configure_for_volume(volume_master_mix)
    p50_single.pick_up_tip()
    for i in range(num_rxns):
        N = len(combinations_to_make[i]['parts'])    # number of part in the combination
        p50_single.consolidate(
            [volume_waterbuffer_per_reaction],          # volume of water and buffer
            [trough.wells_by_name()[well_name] for well_name in ['A1']],
            reaction_plate.wells()[i].bottom(z=0.5), new_tip='never')
    p50_single.drop_tip()

    # This section of the code combines and mix the DNA parts according to the combination list
    # Take one input part and add it in all the combinations it is part of, then do the next input part ...
    p50_single.configure_for_volume(volume_inputDNA)
    for part, combinations in combinations_by_part.items():
        part_well = find_dna(part, dna_plate_map_dict, dna_plate_dict)
        combination_wells = [find_combination(x, combinations_to_make) for x in combinations]
        while combination_wells:
            if len(combination_wells) > 10:
                current_wells = combination_wells[0:10]
                combination_wells = combination_wells[10:]
            else:
                current_wells = combination_wells
                combination_wells = []
            #p50_single.aspirate(volume_inputDNA * len(current_wells), dna_plate_dict[part_well[0]].wells_by_name()[part_well[1]])
            #p50_single.aspirate(volume_inputDNA * len(current_wells), part_well)
            for i in current_wells:
                p50_single.pick_up_tip()
                p50_single.aspirate(volume_inputDNA, part_well)
                p50_single.dispense(volume_inputDNA, i.bottom(z=0.5))
                p50_single.drop_tip()
            #if combination_wells:
                # One washing steps are added to allow recycling of the tips #
                #p50_single.mix(2, 15, water.bottom(z=2))
                #p50_single.blow_out()
        #p50_single.drop_tip()

    protocol.pause('Add the enzyme tube in position B1.')

    p50_single.configure_for_volume(10)
    for i in range(num_rxns):
        p50_single.pick_up_tip()

        p50_single.aspirate(volume_enzyme, well_enzyme.bottom(z=1.5))
        p50_single.dispense(volume_enzyme,  reaction_plate.wells()[i].bottom(z=1))

        p50_single.mix(3, 7, reaction_plate.wells()[i].bottom(z=1))
        p50_single.blow_out()
        p50_single.drop_tip()

    # Seal the Reaction Plate with adhesive film and conduct the GG program
    temp_mod.deactivate()
    temp_mod_reaction.deactivate()
    protocol.pause( 'Please seal the PCR plates and resume run to conduct GG program.')
    '''

    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)
    profile1 = [
            {'temperature': 37, 'hold_time_seconds': 120},
            {'temperature': 16, 'hold_time_seconds': 300}]
    tc_mod.execute_profile(steps=profile1, repetitions=25, block_max_volume=20)
    profile2 = [
            {'temperature': 60, 'hold_time_seconds': 300}]
    tc_mod.execute_profile(steps=profile2, repetitions=1, block_max_volume=20)
    tc_mod.set_block_temperature(4)
    tc_mod.open_lid()
    #temp_mod.set_temperature(4) #Optional
    protocol.pause('Place remove the seal film of the PCR plates and resume run to conduct heat shock program.')
    '''
    temp_mod_reaction.set_temperature(celsius=temp_reaction)
    temp_mod.set_temperature(celsius=temp_reagent)
    protocol.pause('Temperature modules ready, add the competent cells in D1')
    p50_single.configure_for_volume(volume_competent_cells)
    # Add competent cells
    for i in range(0, num_rxns):
            p50_single.pick_up_tip()

            p50_single.aspirate(volume_competent_cells, competent_cell.bottom(z=2), rate =0.2)
            p50_single.dispense(volume_competent_cells, reaction_plate.wells()[i].bottom(z=2), rate =0.2)

            p50_single.mix(1, 25, reaction_plate.wells()[i].bottom(z=2), rate =0.2)
            p50_single.blow_out()
            p50_single.drop_tip()
    temp_mod.deactivate()
    temp_mod_reaction.deactivate()
    protocol.pause('Place seal the PCR plates again and resume run to conduct the heat-shock program.')
    '''
     # Incubate at 4c, then heat shock.
    tc_mod.close_lid()
    profile1 = [
            {'temperature': 4, 'hold_time_seconds': 600},
            {'temperature': 42, 'hold_time_seconds': 90},
            {'temperature': 4, 'hold_time_seconds': 120},
            {'temperature': 37, 'hold_time_seconds': 3600}]
    tc_mod.execute_profile(steps=profile1, repetitions=1, block_max_volume=40)
    tc_mod.set_block_temperature(37)
    tc_mod.open_lid()
    protocol.pause('Please remove the seal and resume for plating')
    '''
    # Plating
    a = len(agar_plate.wells())
    p50_single.configure_for_volume(volume_competent_cells)
    for i in range(0, num_rxns):
        if i <a:
            p50_single.pick_up_tip()
            p50_single.mix(3, volume_competent_cells, reaction_plate.wells()[i].bottom(z=2))
            p50_single.distribute(2.5, reaction_plate.wells()[i].bottom(z=2),
                                  [agar_plate.wells()[i].bottom(z=0).move(position) for position in
                                   [types.Point(x=0, y=0, z=6), types.Point(x=0, y=6, z=5), types.Point(x=6, y=0, z=6),
                                    types.Point(x=0, y=-6, z=5), types.Point(x=-6, y=0, z=6),
                                    types.Point(x=0, y=12, z=5), types.Point(x=7.5, y=7.5, z=6), types.Point(x=12, y=0, z=5),
                                    types.Point(x=7.5, y=-7.5, z=6), types.Point(x=0, y=-12, z=5),
                                    types.Point(x=-7.5, y=-7.5, z=6), types.Point(x=-12, y=0, z=5), types.Point(x=-7.5, y=7.5, z=6)]],
                               disposal_volume=1.5, new_tip='never')
            p50_single.blow_out(trash)
            p50_single.drop_tip()
        if i >= a:
            if i % a == 0:
               protocol.pause('Please change a new agar plates')    
            p50_single.pick_up_tip()
            p50_single.mix(2, volume_competent_cells, reaction_plate.wells()[i].bottom(z=2))
            p50_single.distribute(2.5, reaction_plate.wells()[i].bottom(z=2),
                                  [agar_plate.wells()[i].bottom(z=0).move(position) for position in
                                   [types.Point(x=0, y=0, z=6), types.Point(x=0, y=6, z=5), types.Point(x=6, y=0, z=6),
                                    types.Point(x=0, y=-6, z=5), types.Point(x=-6, y=0, z=6),
                                    types.Point(x=0, y=12, z=5), types.Point(x=7.5, y=7.5, z=6), types.Point(x=12, y=0, z=5),
                                    types.Point(x=7.5, y=-7.5, z=6), types.Point(x=0, y=-12, z=5),
                                    types.Point(x=-7.5, y=-7.5, z=6), types.Point(x=-12, y=0, z=5), types.Point(x=-7.5, y=7.5, z=6)]],
                               disposal_volume=1.5, new_tip='never')
            p50_single.blow_out()
            p50_single.drop_tip()
    #tc_mod.deactivate()
