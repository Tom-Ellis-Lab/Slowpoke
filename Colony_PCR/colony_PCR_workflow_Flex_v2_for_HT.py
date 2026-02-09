
# Colony PCR Protocol V3.0
# Written by Fankang Meng, Imperial College London
#2022-09-04
# Modified by Alicia Da Silva and Henri Galez, Inria and Institut Pasteur
# 2025-08-13

from opentrons import protocol_api, types
import math

metadata = {
    'protocolName': 'Colony PCR - Flex',
    'description': 'Colony PCR'}

requirements = {"robotType": "Flex", "apiLevel": "2.21"}


#####################################
###### PCR reaction settings ########
reaction_volume = 15
dna_volume = 2
enzyme_buffer_volume = 7.5
primer_volume = 1.5
water_volume = 2.5
temperature_modules = 4
#####################################

#number of reactions
num_rxns = len(pcr_recipe_to_make)

def run(protocol: protocol_api.ProtocolContext):
    # Trash need to be specified with Flex
    trash = protocol.load_trash_bin("A3")

    # loading pipette and tips
    tr_50_1 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'C3')
    tr_50_2 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'B3')
    tr_50_3 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'A2')

    p50_single = protocol.load_instrument('flex_1channel_50', 'right', tip_racks=[tr_50_1,tr_50_2,tr_50_3])

    # loading thermocycler
    reaction_mod = protocol.load_module('temperature module gen2', 'A1')
    temp_reaction = reaction_mod.load_adapter('opentrons_96_well_aluminum_block')
    reaction_plate = temp_reaction.load_labware('biorad_96_wellplate_200ul_pcr')

    addition_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'D1')

    # loading plate with picked colonies in 80ul medium
    colony_template_deck= protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'B1')

    # loading rack with PCR recipe tubes
    pcr_mod = protocol.load_module('temperature module gen2', 'D3')
    pcr_deck = pcr_mod.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap')


    pcr_mix_deck = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 'D2')

    # Mettre les modules de température à 4°C
    pcr_mod.set_temperature(temperature_modules)
    reaction_mod.set_temperature(temperature_modules)

    protocol.pause(f'Temp modules ready!')

#Calculate how many PCR reaction systems there are in total
    combinations = [] # list of dict [{name:[...],parts:[water,mastermix,primerfor,primerrev]},{name:[...],parts:[water,mastermix,primerfor,primerrev]}]
    name_1 = pcr_recipe_to_make[0]["name"]
    part_1 = pcr_recipe_to_make[0]["parts"][0:-1]
    combinations.append({"name": [name_1], "parts": part_1})
    for i in pcr_recipe_to_make[1:]:
        name = i["name"]
        part = i["parts"][0:-1] # take water, mastermix, primers (but not colony ofc)
        j = 0
        a = len(combinations) # number of combinations
        while j < a: # looping through existing combinations
            # If the parts of j are different from i["parts"], set a new dictionary in combinations
            name_j = combinations[j]["name"]
            parts_j = combinations[j]["parts"]
            if part == parts_j:
                name_j.append(name)
                break
                # If the parts of j are the same as i["parts"], insert the name corresponding to i["parts"] into the 
                #the dictionary corresponding to this part.
            elif j == a - 1:
                combinations.append({"name": [name], "parts": part})
            j = j + 1

    # This function checks the existance of pcr raw materials and returns for well location of the raw materials
    def find_rawpcr(name, pcr_plate_map_dict, pcr_deck):
        """Return a well containing the named DNA."""
        for plate_name, plate_map in pcr_plate_map_dict.items():
            for i, row in enumerate(plate_map):
                for j, dna_name in enumerate(row):
                    if dna_name == name:
                        well_num = 4 * j + i
                        return pcr_deck.wells()[well_num]
        raise ValueError("Could not find dna piece named \"{0}\"".format(name))

    # This function checks if the DNA parts exist in the DNA plates and returns for well location of output DNA combinations
    def find_combination(name, combinations_to_make):
        """Return a well containing the named combination."""
        for i, combination in enumerate(combinations_to_make):
            if combination["name"] == name:
                if i < 96:
                    return reaction_plate.wells()[i]
                else:
                    return addition_plate.wells()[i-96]
        raise ValueError("Could not find combination \"{0}\".".format(name))

    #According to the type of PCR reaction, add different PCR raw materials and distribute them into the corresponding locations.
    for i, combination in enumerate(combinations):
        name_i = combination["name"]
        part_i = combination["parts"]
        pcr_sample_number = len(name_i) * 1.2 # make for 20% extra samples to avoid pipetting error
        pcr_plate_map_dict_list = pcr_deck_colony_template_maps_dict["pcr_deck_map"]
        pcr_plate_map_dict_number = sum([len(count) for count in pcr_plate_map_dict_list])

        for j, part in enumerate(part_i):
            if j == 0:
                a = water_volume
            elif j == 1:
                a = enzyme_buffer_volume
            elif j == 2 or 3:
                a = primer_volume

            volume_j = pcr_sample_number * a
            repeat = volume_j // 50
            last = volume_j % 50

            rawpcr_well = find_rawpcr(part, pcr_deck_colony_template_maps_dict, pcr_deck)

            p50_single.pick_up_tip()
            for k in range(int(repeat)):
                p50_single.configure_for_volume(50)
                p50_single.transfer(50,
                                     rawpcr_well.bottom(z=1),
                                     pcr_mix_deck.wells()[i].bottom(z=2),
                                     blow_out=True, blowout_location='destination well',
                                     new_tip='never')
                if (j == 2) or (j == 3):
                    p50_single.drop_tip()
                    p50_single.pick_up_tip()
                elif (k % 4) == 3:
                    p50_single.drop_tip()
                    p50_single.pick_up_tip()

            p50_single.configure_for_volume(last)
            p50_single.transfer(last,
                                 rawpcr_well.bottom(z=1),
                                 pcr_mix_deck.wells()[i].bottom(z=2),
                                 blow_out=True, blowout_location='destination well',
                                 new_tip='never')
            p50_single.drop_tip()

        protocol.pause('Mix PCR mastermixes manually if needed')

        p50_single.pick_up_tip()
        volume_mix = min(50,(reaction_volume - dna_volume) * (pcr_sample_number - 1))
        p50_single.configure_for_volume(volume_mix)
        p50_single.mix(2, volume_mix, pcr_mix_deck.wells()[i].bottom(z=1))
        p50_single.drop_tip()

        p50_single.configure_for_volume(reaction_volume-dna_volume)
        pcr_combination_wells = [find_combination(x, pcr_recipe_to_make) for x in name_i]

        nb_per_disp = 3 * (50 // reaction_volume-dna_volume)
        div = len(pcr_combination_wells) // nb_per_disp
        for disp in range(div + 1):
            start_pos = disp * nb_per_disp
            end_pos = min(start_pos + nb_per_disp, len(pcr_combination_wells))
            distribute_wells = pcr_combination_wells[start_pos:end_pos]
            if distribute_wells != []:
                p50_single.distribute(reaction_volume-dna_volume,
                                pcr_mix_deck.wells()[i].bottom(z=1),
                                distribute_wells,
                                disposal_volume=1, new_tip='once')


    # This function checks the existence of pcr raw materials and returns for well location of the raw materials
    def find_template(name, pcr_deck_colony_template_maps_dict, colony_template_deck):
        """Return a well containing the named DNA."""
        for plate_name, plate_map in pcr_deck_colony_template_maps_dict.items():
            for i, row in enumerate(plate_map):
                for j, colony_name in enumerate(row):
                    if colony_name == name:
                        well_num = 8 * j + i
                        return colony_template_deck.wells()[well_num]
        raise ValueError("Could not find colony template named \"{0}\"".format(name))

    combinations_by_colony_template = {}
    for i in pcr_recipe_to_make:
        name = i["name"]
        template = i["parts"][-1]
        if template in combinations_by_colony_template.keys():
            combinations_by_colony_template[template].append(name)
        else:
            combinations_by_colony_template[template] = [name]

    p50_single.configure_for_volume(dna_volume)
    for part, combination_template in combinations_by_colony_template.items():
        template_well = find_template(part, pcr_deck_colony_template_maps_dict, colony_template_deck)
        colony_combination_wells = [find_combination(x, pcr_recipe_to_make) for x in combination_template]
        for colony_well in colony_combination_wells:
            p50_single.pick_up_tip()
            p50_single.aspirate(dna_volume, template_well.bottom(z=2))
            p50_single.dispense(dna_volume, colony_well)
            mix_volume = min(reaction_volume * 0.75, 10)
            p50_single.mix(3, mix_volume, colony_well)
            p50_single.blow_out()
            p50_single.drop_tip()



    # seal the pcr plate with adhesive film and conduct the GG program
    protocol.pause('Please seal the PCR plates.')
    pcr_mod.deactivate()
    reaction_mod.deactivate()

