# Colony PCR Protocol V3.0
# Written by Fankang Meng and Koray Malci Imperial College London
#2022-09-04
# Modified by Alicia Da Silva and Henri Galez, Inria and Institut Pasteur
# 2025-08-13

from opentrons import protocol_api, types
import math

metadata = {
    'protocolName': 'Colony_PCR_Flex',
    'description': 'Colony PCR'}

requirements = {"robotType": "Flex", "apiLevel": "2.21"}


#####################################
#PCR cycles and extension time；change this!!!
pcr_cycles = [
    {'temperature': 95, 'hold_time_seconds': 30},
    {'temperature': 55 , 'hold_time_seconds': 30},
    {'temperature': 72, 'hold_time_seconds': 60}]#change this!!!
pcr_cycles_number = 30 #change this!!!
reaction_volume = 10 #change this!!!
#####################################


#number of ractions
num_rxns = len(pcr_recipe_to_make)

def run(protocol: protocol_api.ProtocolContext):
    # Trash need to be specified with Flex
    trash = protocol.load_trash_bin("A3")

    # loading pipette and tips
    #tr_1000 = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'B3')
    tr_50 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', 'C3')

    p50_single = protocol.load_instrument('flex_1channel_50', 'right', tip_racks=[tr_50])
    #p1000_single = protocol.load_instrument('flex_1channel_1000', 'left', tip_racks=[tr_1000])

    #p10_single = protocol.load_instrument('p10_single', 'right', tip_racks=[tr_20])
    #p300_single = protocol.load_instrument('p300_single', 'left', tip_racks=[tr_300])

    # loading temperature module
    # tc_mod = protocol.load_module('Thermocycler Module')
    # reaction_plate = tc_mod.load_labware('biorad_96_wellplate_200ul_pcr')
    #temp_mod_pcr = protocol.load_module('temperature module gen2', 'A1')
    #reaction_plate = temp_mod_pcr.load_labware('biorad_96_wellplate_200ul_pcr', 'A1')
    reaction_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'C1', 'colony_template_plate')

    addition_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 'D1')
    
    # loading plate with picked colonies in 80ul medium
    colony_template_deck= protocol.load_labware('biorad_96_wellplate_200ul_pcr', 'D2', 'colony_template_plate')

    # loading rack with PCR recipe tubes
    pcr_mod = protocol.load_module('temperature module gen2', 'D3')
    pcr_deck = pcr_mod.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap', 'D3')
    
    pcr_mix_mod = protocol.load_module('temperature module gen2', 'A1')
    pcr_mix_deck = pcr_mix_mod.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap', 'A1')
    
    # Mettre les modules de température à 8°C
    pcr_mod.set_temperature(8)
    pcr_mix_mod.set_temperature(8)
    print("Modules de temperature regles a 8°C")

    # Optionnel : attendre que la température soit atteinte
    protocol.comment("En attente que les modules atteignent 8°C...")
    pcr_mod.await_temperature(8)
    pcr_mix_mod.await_temperature(8)
    protocol.comment("Temperature de 8°C atteinte")

    #pcr_deck = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '9', ' pcr_deck')
    #pcr_mix_deck = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '5', ' pcr_deck')

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

    # This function checks if the DNA parts exist in the DNA plates and returns for well locaion of output DNA combinations
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
    
    if len(combinations) > 0:
        combination = combinations[0]
        part_i = combination["parts"]
        total_reactions = len(pcr_recipe_to_make)
    
        #Prepare Dreamtaq master mix for all reactions
        for j, part in enumerate(part_i):
            if j == 0:
                a = 4.5  # Water - 4.5 µL per reaction
            elif j == 1:
                a = 2.5  # MM_Dreamtaq - 2.5 µL per reaction
            elif j == 2:
                a = 1    # P128 - 1 µL per reaction
            elif j == 3:
                a = 1    # P132 - 1 µL per reaction
            
            rawpcr_well = find_rawpcr(part, pcr_deck_colony_template_maps_dict, pcr_deck)
            p50_single.pick_up_tip()
            p50_single.transfer(total_reactions * a + 2,  # +2 of dead volume
                            rawpcr_well.bottom(z=2),
                            pcr_mix_deck.wells()[0].bottom(z=3),
                            blow_out=True, blowout_location='destination well',
                            new_tip='never')
            p50_single.drop_tip()

        # Mix the master mix
        p50_single.pick_up_tip()
        p50_single.mix(3, 40, pcr_mix_deck.wells()[0].bottom(z=1))
        p50_single.drop_tip()

        # Distribute the master mix
        destination_wells = []
        for idx in range(total_reactions):
            if idx < 96:
                destination_wells.append(reaction_plate.wells()[idx])
            else:
                destination_wells.append(addition_plate.wells()[idx-96])

        p50_single.pick_up_tip()
        p50_single.distribute(reaction_volume-1,
                            pcr_mix_deck.wells()[0].bottom(z=1),
                            destination_wells,
                            disposal_volume=2, 
                            new_tip='never')
        p50_single.drop_tip()

        print(f"Master mix distribué dans {total_reactions} puits")

    # This function checks the existance of pcr raw materials and returns for well location of the raw materials
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

    # Transfert the colonies in the different wells 
    for idx, recipe in enumerate(pcr_recipe_to_make):
        colony_name = recipe["parts"][-1]  # The colonie is the last element 
                
        template_well = find_template(colony_name, pcr_deck_colony_template_maps_dict, colony_template_deck)
        
        if idx < 96:
            destination_well = reaction_plate.wells()[idx]
        else:
            destination_well = addition_plate.wells()[idx-96]
    
        p50_single.pick_up_tip()
        p50_single.transfer(1,
                            template_well.bottom(z=1),
                            destination_well.bottom(z=1),
                            blow_out=True, 
                            blowout_location='destination well',
                            new_tip='never')
        p50_single.drop_tip()
        
        print(f"Transfere colonie {colony_name} vers puits {idx} (plaque {'reaction' if idx < 96 else 'addition'})")
        
        
    # Turn off the modules
    pcr_mod.deactivate()
    pcr_mix_mod.deactivate()
    print("Modules de temperature desactives")

    '''

    # seal the pcr plate with adhesive film and conduct the GG program
    protocol.pause( 'Please seal the PCR plates and resume run to conduct PCR program.')
    
    # PCR steps
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)
    # denaturation
    pre_denaturation_step = [{'temperature': 98, 'hold_time_seconds': 300}]
    tc_mod.execute_profile(steps=pre_denaturation_step, repetitions=1, block_max_volume=reaction_volume)
    # amplification cycles
    j = 0
    while j < pcr_cycles_number:
        tc_mod.execute_profile(steps=pcr_cycles, repetitions=1, block_max_volume=reaction_volume)
        j = j + 1
    # further  extension step
    extension_step = [{'temperature': 72, 'hold_time_seconds': 300}]
    tc_mod.execute_profile(steps=extension_step, repetitions=1, block_max_volume=reaction_volume)
    tc_mod.set_block_temperature(25)


    tc_mod.open_lid()
    tc_mod.deactivate()
    '''
