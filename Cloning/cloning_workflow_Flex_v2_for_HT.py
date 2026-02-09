

# YTK/STK/KTK GG, Transformation, and Plating protocol
# Written by Fankang Meng, Imperial College London
# Updated by Alicia Da Silva and Henri Galez for Flex robot, Institut Pasteur

from opentrons import protocol_api, types
import time
import math


metadata = {
    'protocolName': 'Golden Gate Cloning - Flex',
    'description': 'GG & Transformation & plating using a Flex robot for genetic toolkits.'}

requirements = {"robotType": "Flex", "apiLevel": "2.21"}

num_rxns = len(combinations_to_make)
volume_buffer = 1.2
volume_enzyme = 1.2
volume_reaction = 12
nb_parts = 6
volume_inputDNA = 1
volume_waterbuffer_per_reaction = volume_reaction - volume_buffer - nb_parts * volume_inputDNA
volume_master_mix = volume_enzyme + volume_waterbuffer_per_reaction

volume_tubes_competent = 1100
volume_tubes_competent_safe = volume_tubes_competent - 100
volume_competent_cells = 50
nb_reaction_per_tube = volume_tubes_competent_safe // volume_competent_cells

temp_reaction = 4
temp_reagent = 4

def run(protocol: protocol_api.ProtocolContext):

    # Compute needed tips
    def calculate_tips_needed():
        # 1. Tips for buffer/water: 1 tip total
        nb_per_disp = 2 * (50 // volume_waterbuffer_per_reaction)  # number of wells that can be distributed per dispense (2 distribute per tip)
        buffer_tips = num_rxns // nb_per_disp

        # 2. Tips for DNA parts: 1 tip per parts per combination
        dna_tips = 0
        combinations_by_part = {}
        for combo in combinations_to_make:
            for part in combo["parts"]:
                if part in combinations_by_part:
                    combinations_by_part[part] += 1
                else:
                    combinations_by_part[part] = 1
        dna_tips = sum(combinations_by_part.values())

        # 3. Tips for enzyme: 1 tip per reaction
        enzyme_tips = num_rxns

        # 4. Tips for competent cells: 1 tip per reaction
        competent_tips = num_rxns

        # 5. Tips for plating: 1 tip per reaction
        plating_tips = num_rxns

        total_tips = buffer_tips + dna_tips + enzyme_tips + competent_tips + plating_tips

        # Add a 10% safety margin
        total_tips = int(total_tips * 1.1)

        return total_tips, {
            'buffer': buffer_tips,
            'dna': dna_tips,
            'enzyme': enzyme_tips,
            'competent': competent_tips,
            'plating': plating_tips
        }

    # Calculation of reagent quantities
    total_buffer_needed = volume_waterbuffer_per_reaction * num_rxns
    total_enzyme_needed = volume_enzyme * num_rxns
    total_competent_cells_needed = volume_competent_cells * num_rxns

    tips_needed, tips_breakdown = calculate_tips_needed()
    tips_per_rack = 96
    racks_needed = math.ceil(tips_needed / tips_per_rack)

    # Slots available for tip racks
    available_slots = ['A2', 'B1', 'B2', 'D1','C3','D2']  # Emplacements libres

    # Pause for tip rack setup
    setup_message = f""" Tip setup:
- Number of constructions: {num_rxns}
- Tips needed : {tips_needed}
- Tip racks needed: {racks_needed}

Place {racks_needed} of 50 uL at the location :
"""

    for i in range(racks_needed):
        if i < len(available_slots):
            setup_message += f"\n - Rack {i+1} of 50uL: {available_slots[i]}"

    protocol.pause(setup_message)

    # Trash need to be specified with Flex
    trash = protocol.load_trash_bin("A3")

    # Dynamic tip rack loading
    tip_racks = []

    # Existing racks
    tr_300 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B3', 'Tips Rack 1')
    tip_racks.extend([tr_300])

    # Additional racks if needed
    for i in range(min(racks_needed - 1, len(available_slots))):  # -1 because we already have 1 tip rack
        if racks_needed > 1:  # Only if we need more than one rack
            additional_rack = protocol.load_labware('opentrons_flex_96_tiprack_50ul',
                                                   available_slots[i],
                                                   f'Tips Rack {i+3}')
            tip_racks.append(additional_rack)

    # Load in pipettes
    p50_single = protocol.load_instrument('flex_1channel_50', 'right', tip_racks=tip_racks)

    # Load modules
    temp_mod_reaction = protocol.load_module('temperature module gen2', 'A1')
    temp_mod_reaction.set_temperature(celsius=temp_reaction)
    temp_adapter = temp_mod_reaction.load_adapter('opentrons_96_well_aluminum_block')
    reaction_plate = temp_adapter.load_labware('biorad_96_wellplate_200ul_pcr')

    temp_mod = protocol.load_module('temperature module gen2', 'D3')
    temp_mod.set_temperature(celsius=temp_reagent)
    trough = temp_mod.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap', 'D3')
    buffer_H2O_mix = trough.wells()[0]  # Well A1
    well_enzyme = trough.wells()[1]  # Well B1
    dilution_water = trough.wells()[2]  # Well C1
    #competent_cell = trough.wells()[3]  # Well D1
    competent_cells = [trough.wells()[3], trough.wells()[7], trough.wells()[11], trough.wells()[15], trough.wells()[19]]  # Well D1 -> D5
    liquid_waste = trough.wells()[4]  # Well A2

    # Load in Input DNA Plate

    dna_plate_dict = {}
    plate_name = list(dna_plate_map_dict.keys())
    dna_plate_dict[plate_name[0]] = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'C2', 'Input DNA Plate')
    #activate the following line by deleting '#' if a second custom_parts_map needs to be used
    #dna_plate_dict[plate_name[1]] = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'D2', 'Input DNA Plate2')

    # Load in Agar plate
    agar_plate = protocol.load_labware('corning_6_wellplate_16.8ml_flat', 'C1', 'Agar Plate')


    # This function checks the existance of DNA parts and returns for well location of the parts
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
                        if 'fixed_toolkit_map' in plate_name:
                            #well_num = 8 * j + i
                            well_num = rows[i]+str(columns[j])
                        elif 'custom_parts_map' in plate_name:
                            #well_num = 4 * j + i
                            well_num = rows[i]+str(columns[j])
        if well_num is not None:
            return dna_plate_dict[plate_part].wells_by_name()[well_num]
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

    # Step 1: Add Buffer/Water
    protocol.pause(f'Temperature modules ready!\n Put {total_buffer_needed * 1.2} uL of buffer/water in A1 position.')

    p50_single.configure_for_volume(volume_waterbuffer_per_reaction)

    nb_per_disp = 2 * (30 // math.ceil(volume_waterbuffer_per_reaction))  # number of wells that can be distributed per dispense (2 distribute per tip)
    div = num_rxns // nb_per_disp
    for disp in range(div + 1):
        start_pos = disp * nb_per_disp
        end_pos = min(start_pos + nb_per_disp, num_rxns)
        distribute_wells = reaction_plate.wells()[start_pos:end_pos]
        if distribute_wells !=[]:
            p50_single.distribute(volume_waterbuffer_per_reaction,
                                  [trough.wells_by_name()[well_name] for well_name in ['A1']],
                                  distribute_wells,
                                  disposal_volume=1, new_tip='always')

    # Step 2: Add DNA parts
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
            for i in current_wells:
                p50_single.pick_up_tip()
                p50_single.aspirate(volume_inputDNA, part_well.bottom(z=1))
                p50_single.dispense(volume_inputDNA, i.bottom(z=1))
                p50_single.drop_tip()

    # Step 3: Add enzyme
    protocol.pause(f'Put {total_enzyme_needed} uL of enzyme in B1')

    p50_single.configure_for_volume(10)
    for i in range(num_rxns):
        p50_single.pick_up_tip()
        p50_single.aspirate(volume_enzyme, well_enzyme.bottom(z=1.5))
        p50_single.dispense(volume_enzyme,  reaction_plate.wells()[i].bottom(z=1))
        mix_volume = min(volume_reaction*0.75, 10)
        p50_single.mix(3, mix_volume, reaction_plate.wells()[i].bottom(z=1))
        p50_single.blow_out()
        p50_single.drop_tip()

    # Step 4 : Incubation GG
    temp_mod.deactivate()
    temp_mod_reaction.deactivate()
    protocol.pause('Golden Gate:\nn Seal PCR plates with adhesive film\n Start the Golden Gae program (cycles 37C/16C)\n Press Resume once finished.')

    temp_mod_reaction.set_temperature(celsius=temp_reaction)
    temp_mod.set_temperature(celsius=temp_reagent)

    # Step 5: Add competent cells
    protocol.pause(f'{total_competent_cells_needed} uL of competent cells in total, {volume_tubes_competent} per tube in D1 -> D5')

    p50_single.configure_for_volume(volume_competent_cells)
    for i in range(0, num_rxns):
        tube_number = i // nb_reaction_per_tube
        competent_cell = competent_cells[tube_number]
        p50_single.pick_up_tip()
        p50_single.aspirate(volume_competent_cells, competent_cell.bottom(z=2), rate =0.2)
        p50_single.dispense(volume_competent_cells, reaction_plate.wells()[i].bottom(z=2), rate =0.2)
        p50_single.mix(1, 25, reaction_plate.wells()[i].bottom(z=2), rate =0.2)
        p50_single.blow_out()
        p50_single.drop_tip()

    temp_mod.deactivate()
    temp_mod_reaction.deactivate()

    # Step 6: heat shock
    protocol.pause(' Heat shock:\n Reseal the PCR plates\n Proceed with the heat shock program \n Press Resume to begin plating.')

    # Step 7: plating
    num_agar_plates_needed = math.ceil(num_rxns / 6)
    total_volume_per_construct = 2.5 * 13  #13 deposition points per construct
    total_plating_volume = total_volume_per_construct * num_rxns

    plating_setup_message = f""" Setup plating:
 {num_rxns} constructions to plate
 {num_agar_plates_needed} agar plaque(s)
 Total volume to plate: {total_plating_volume} uL

Place the first agar plate in position C1 and press Resume."""

    protocol.pause(plating_setup_message)

    p50_single.configure_for_volume(volume_competent_cells)
    wells_per_plate = 6

    for i in range(0, num_rxns):
        well_index = i % wells_per_plate

        if well_index == 0 and i > 0:
            plate_number = (i // wells_per_plate) + 1
            protocol.pause(f' Changing agar plate:\n Remove the full agar plate (plate {plate_number - 1})\n Place a new empty agar plate at the same location C1\n You start the plate {plate_number}/{num_agar_plates_needed}\nPress Resume once the new plate is in place.')

        current_well = agar_plate.wells()[well_index]

        positions = [
            types.Point(x=0, y=0, z=6), types.Point(x=0, y=6, z=5), types.Point(x=6, y=0, z=6),
            types.Point(x=0, y=-6, z=5), types.Point(x=-6, y=0, z=6),
            types.Point(x=0, y=12, z=5), types.Point(x=7.5, y=7.5, z=6), types.Point(x=12, y=0, z=5),
            types.Point(x=7.5, y=-7.5, z=6), types.Point(x=0, y=-12, z=5),
            types.Point(x=-7.5, y=-7.5, z=6), types.Point(x=-12, y=0, z=5), types.Point(x=-7.5, y=7.5, z=6)
        ]

        p50_single.pick_up_tip()
        p50_single.mix(3, volume_competent_cells, reaction_plate.wells()[i].bottom(z=2))
        p50_single.distribute(2.5, reaction_plate.wells()[i].bottom(z=2),
                            [current_well.bottom(z=0).move(position) for position in positions],
                            disposal_volume=1.5, new_tip='never')
        p50_single.blow_out(trash)
        p50_single.drop_tip()

    # Final message
    final_message = f""" PROTOCOL COMPLETED!

 NEXT STEPS:
- Remove the last agar plate
- Incubate the agar plates at 37C overnight
- Check colony growth tomorrow

Congrats! """

    protocol.pause(final_message)