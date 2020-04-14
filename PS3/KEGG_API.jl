function download_reactions_from_kegg(path_to_ec_array,path_to_reaction_file)

    # initalize -
    record_array = String[]

    # load the ec_array file -
    ec_array = readdlm(path_to_ec_array,',')

    # how many ec's do we have?
    (number_of_ec_numbers,nc) = size(ec_array)

    for ec_number_index = 1:number_of_ec_numbers

        # get the ec_number -
        rxn_index = ec_array[ec_number_index,1]
        ec_number = ec_array[ec_number_index,2]

        # use the KEGG API to download the reaction -
        ec_rn_linkage = readstring(`curl -X GET http://rest.kegg.jp/link/rn/$(ec_number)`)

        # if the server returned something, it should be in the form ec# rn#
        # we need the rn# to pull down the reaction string -
        if (length(ec_rn_linkage)>1)

            # ec_rn_linkage can have more than one association -
            # split around the \n, and then process each item
            P = reverse(split(chomp(ec_rn_linkage),'\n'))
            while (!isempty(P))

                # pop -
                local_record = pop!(P)



                # The second fragment is the rn# -
                # split around the \t
                rn_number = string(split(local_record,'\t')[2])

                # ok, run a second query to pull the reaction string down -
                reaction_string = chomp(readstring(`curl -X GET http://rest.kegg.jp/list/$(rn_number)`))

                # if we have something in reaction string, create a string record -
                if (length(reaction_string)!=0)
                    tmp_string = "$(rxn_index),$(ec_number),$(reaction_string)"
                    push!(record_array,tmp_string)
                end
            end
        end
    end

    # dump reaction array to disk -
    open(path_to_reaction_file, "w") do f

        for line_item in record_array
            write(f,"$(line_item)\n")
        end
    end
end
