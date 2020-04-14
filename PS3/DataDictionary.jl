function check_unbalanced_boundary_reactions()

	# load the original data dictionary -
	data_dictionary = maximize_urea_production()

	# add boundary species to the list of species -
	list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]
	list_of_metabolite_symbols = [
		list_of_metabolite_symbols ;

		"M_Carbamoyl_phosphate_b"	;	# 1
		"M_L-Aspartate_b"	;			# 2
		"M_Fumarate_b"	;				# 3
		"M_Urea_b"	;					# 4
		"M_ATP_b"	;					# 5
		"M_AMP_b"	;					# 6
		"M_Diphosphate_b"	;			# 7
		"M_Orthophosphate_b"	;		# 8
		"M_Oxygen_b"	;				# 9
		"M_NADPH_b"	;					# 10
		"M_H_b"	;						# 11
		"M_Nitric_oxide_b"	;			# 12
		"M_NADP_b"	;					# 13
		"M_H2O_b"	;					# 14
	]

	# add some boundary species (rows to the stoichiometric array)
	# there will be in the same order as the list_of_metabolite_symbols
	S = data_dictionary["stoichiometric_matrix"]

	# transfer reaction block -
	transfer_reaction_block = zeros(14,21)
	transfer_reaction_block[1,7] = -1 		# 1 M_Carbamoyl_phosphate_b -> box
	transfer_reaction_block[2,8] = -1 		# 2 M_L-Aspartate_b -> box
	transfer_reaction_block[3,9] =	1		# 3 box -> M_Fumarate_b
	transfer_reaction_block[4,10] = 1		# 4 box -> M_Urea_b
	transfer_reaction_block[5,11] = -1		# 5 M_ATP_b -> box
	transfer_reaction_block[6,12] = 1		# 6 box -> M_AMP_b
	transfer_reaction_block[7,13] = 1		# 7 box -> M_Diphosphate_b
	transfer_reaction_block[8,14] = 1		# 8 box -> M_Orthophosphate_b
	transfer_reaction_block[9,15] = -1		# 9 M_Oxygen_b -> box
	transfer_reaction_block[10,16] = -1		# 10 M_NADPH_b -> box
	transfer_reaction_block[11,17] = -1		# 11 M_H_b -> box
	transfer_reaction_block[12,18] = 1		# 12 box -> M_Nitric_oxide_b
	transfer_reaction_block[13,19] = 1		# 13 box -> M_NADP_b
	transfer_reaction_block[14,20] = 1		# 14 box -> M_H2O_b
	transfer_reaction_block[14,21] = -1		# 15 M_H2O_b -> box
	SAUG = [S ; transfer_reaction_block]

	# repackage -
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["stoichiometric_matrix"] = SAUG

	# return the modified data dictionary -
	return data_dictionary
end



function maximize_urea_production_open()

	# load the original data dictionary -
	data_dictionary = maximize_urea_production()

	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
	default_flux_bounds_array[15,1] = 0
	default_flux_bounds_array[15,2] = 10

	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
	default_flux_bounds_array[16,1] = 0
	default_flux_bounds_array[16,2] = 10

	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
	default_flux_bounds_array[17,1] = 0
	default_flux_bounds_array[17,2] = 10

	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
	default_flux_bounds_array[18,1] = -10
	default_flux_bounds_array[18,2] = 10

	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
	default_flux_bounds_array[19,1] = -10
	default_flux_bounds_array[19,2] = 10

	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
	default_flux_bounds_array[20,1] = -10
	default_flux_bounds_array[20,2] = 10

	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	default_flux_bounds_array[21,1] = 0
	default_flux_bounds_array[21,2] = 10

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

	# return the updated dictionary -
	return data_dictionary
end


function maximize_urea_production()

	# load the original data dictionary -
	data_dictionary = DataDictionary()

	# 1: set the objective function -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[10] = -1

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# let all exchanges go from 0,10 mmol/gDW-hr
	range_of_exchange_reactions = collect(7:21)
	for reaction_index in range_of_exchange_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = 10.0
	end

	# don't allow water exchange -
	default_flux_bounds_array[20,2] = 0
	default_flux_bounds_array[21,2] = 0

	# we have some specific values for v1 -> v5
	E = 0.01E-3	# mmol/gDW
	metabolic_vmax_array = [
		203*(3600)*E	;	# v1 ec:6.3.4.5 mmol/gDW-hr
		34.5*(3600)*E	;	# v2 ec:4.3.2.1 mmol/gDW-hr
		249*(3600)*E	;	# v3 ec:3.5.3.1 mmol/gDW-hr
		88.1*(3600)*E	;	# v4 ec:2.1.3.3 mmol/gDW-hr
		13.7*(3600)*E	;	# v5 ec:1.14.13.39 mmol/gDW-hr
		13.7*(3600)*E	;	# v6 ec:1.14.13.39 mmol/gDW-hr
	]
	range_of_cycle_reactions = collect(1:6)
	for reaction_index in range_of_cycle_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = metabolic_vmax_array[reaction_index]
	end

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array

	# return the updated dictionary -
	return data_dictionary
end


function DataDictionary()

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	# What is the system dimension? -
	(number_of_species,number_of_reactions) = size(stoichiometric_matrix)

	# Setup default flux bounds array -
	default_bounds_array = [
		0	metabolic_vmax_array[1]	;	# Vmax [mmol/gdw-hr] 1	M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
		0	metabolic_vmax_array[2]	;	# Vmax [mmol/gdw-hr] 2	M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
		0	metabolic_vmax_array[3]	;	# Vmax [mmol/gdw-hr] 3	M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
		0	metabolic_vmax_array[4]	;	# Vmax [mmol/gdw-hr] 4	M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
		0	metabolic_vmax_array[5]	;	# Vmax [mmol/gdw-hr] 5	2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
		0	metabolic_vmax_array[6]	;	# Vmax [mmol/gdw-hr] 6	2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c

		0 10	;	# Vmax [mmol/gdw-hr] 7	[] --> M_Carbamoyl_phosphate_c
		0 10	;	# Vmax [mmol/gdw-hr] 8	[] --> M_L-Aspartate_c
		-10 10	;	# Vmax [mmol/gdw-hr] 9	M_Fumarate_c --> []
		-10 10	;	# Vmax [mmol/gdw-hr] 10	M_Urea_c --> []
		0 10	;	# Vmax [mmol/gdw-hr] 11	[] --> M_ATP_c
		-10 10	;	# Vmax [mmol/gdw-hr] 12	M_AMP_c --> []
		-10 10	;	# Vmax [mmol/gdw-hr] 13	M_Diphosphate_c --> []
		-10 10	;	# Vmax [mmol/gdw-hr] 14	M_Orthophosphate_c --> []
		0 10	;	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
		0 10	;	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
		0 10	;	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
		-10 10	;	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
		-10 10	;	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
		-10 10	;	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
		0 10	;	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
	];

	# Setup default species bounds array -
		species_bounds_array = [
			0.0	0.0	;	# 1 M_AMP_c
			0.0	0.0	;	# 2 M_ATP_c
			0.0	0.0	;	# 3 M_Carbamoyl_phosphate_c
			0.0	0.0	;	# 4 M_Diphosphate_c
			0.0	0.0	;	# 5 M_Fumarate_c
			0.0	0.0	;	# 6 M_H2O_c
			0.0	0.0	;	# 7 M_H_c
			0.0	0.0	;	# 8 M_L-Arginine_c
			0.0	0.0	;	# 9 M_L-Aspartate_c
			0.0	0.0	;	# 10 M_L-Citrulline_c
			0.0	0.0	;	# 11 M_L-Ornithine_c
			0.0	0.0	;	# 12 M_N-(L-Arginino)succinate_c
			0.0	0.0	;	# 13 M_NADPH_c
			0.0	0.0	;	# 14 M_NADP_c
			0.0	0.0	;	# 15 M_Nitric_oxide_c
			0.0	0.0	;	# 16 M_Orthophosphate_c
			0.0	0.0	;	# 17 M_Oxygen_c
			0.0	0.0	;	# 18 M_Urea_c
		];

	# Min/Max flag - default is minimum -
		is_minimum_flag = true

		# Metabolic Vmax array (units: mmol/B-hr) -
		v1 = v1_max*((4.67E-3)/(3.92E-4+4.67E-3))*((1.49E-2)/(1.54E-4+1.49E-2)) #6.3.4.5
		v2 = v2_max*(1) #4.3.2.1
		v3 = v3_max*((2.55E-4)/(1.55E-3+2.55E-4)) #3.5.3.1
		v4 = v4_max*((1.01E-5)/(8.50E-4+1.01E-5)) #2.1.3.3
		v5 = v5_max*((2.55E-4)/(3.50E-6+2.55E-4))#1.14.13.39
		v6 = v6_max*(1)
		v_b = 10000/3600 #mmol/gdw-hr

		metabolic_vmax_array = [
			v1	;	# Vmax [mmol/gdw-hr] 1	M_ATP_c+M_L-Citrulline_c+M_L-Aspartate_c --> M_AMP_c+M_Diphosphate_c+M_N-(L-Arginino)succinate_c
			v2	;	# Vmax [mmol/gdw-hr] 2	M_N-(L-Arginino)succinate_c --> M_Fumarate_c+M_L-Arginine_c
			v3	;	# Vmax [mmol/gdw-hr] 3	M_L-Arginine_c+M_H2O_c --> M_L-Ornithine_c+M_Urea_c
			v4	;	# Vmax [mmol/gdw-hr] 4	M_Carbamoyl_phosphate_c+M_L-Ornithine_c --> M_Orthophosphate_c+M_L-Citrulline_c
			v5	;	# Vmax [mmol/gdw-hr] 5	2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c --> 2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c
			v6 	;	# Vmax [mmol/gdw-hr] 6	2.0*M_Nitric_oxide_c+2.0*M_L-Citrulline_c+3.0*M_NADP_c+4.0*M_H2O_c --> 2.0*M_L-Arginine_c+4.0*M_Oxygen_c+3.0*M_NADPH_c+3.0*M_H_c

			v_b	;	# Vmax [mmol/gdw-hr] 7	[] --> M_Carbamoyl_phosphate_c
			v_b	;	# Vmax [mmol/gdw-hr] 8	[] --> M_L-Aspartate_c
			v_b	;	# Vmax [mmol/gdw-hr] 9	M_Fumarate_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 10	M_Urea_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 11	[] --> M_ATP_c
			v_b	;	# Vmax [mmol/gdw-hr] 12	M_AMP_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 13	M_Diphosphate_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 14	M_Orthophosphate_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 15	[] --> M_Oxygen_c
			v_b	;	# Vmax [mmol/gdw-hr] 16	[] --> M_NADPH_c
			v_b	;	# Vmax [mmol/gdw-hr] 17	[] --> M_H_c
			v_b	;	# Vmax [mmol/gdw-hr] 18	M_Nitric_oxide_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 19	M_NADP_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 20	M_H2O_c --> []
			v_b	;	# Vmax [mmol/gdw-hr] 21	[] --> M_H2O_c
		];


			cell_mass = 2.8E-13     # g
			cell_volume = 1.0E-12;  #L
			water_fraction = 0.798;
			cell_drymass = (1 - water_fraction)*cell_mass;

	return data_dictionary
end
