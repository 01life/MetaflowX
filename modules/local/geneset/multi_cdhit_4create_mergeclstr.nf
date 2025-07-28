
process CDHITCLSTR {

    label 'process_single'

    input:
    path(clstr)
    path(div_o)
    path(task_num)


    output:
    path("clstr_order.txt"),emit:"clstr_order"
    path("unique.fa"),emit:"unique_fa"

    script:
    def split_num = task_num.getSimpleName().toInteger()

    """

    # Initialize the reps array.
    reps=()
    seg_no=${split_num}
    out_clstr=clstr_order.txt
    out=unique.fa
    indiv='all.cds.fa.div'

    #Iterate \$i from 0 to \$seg_no.
    for (( i=0; i<seg_no; i++ ))
    do
            master_clstr="\$indiv-\$i-o.clstr"

            # If the master_clstr file does not exist or is empty, skip the current loop.
            [ -s "\$master_clstr" ] || continue

            this_rep="\$indiv-\$i-o"

            # If the this_rep file does not exist, output an error message and exit.
            if [ ! -e "\$this_rep" ]; then
                echo "No rep \$this_rep"
                exit 1
            fi

            # Add this_rep to the reps array.
            reps+=("\$this_rep")

            # Initialize the slave_clstr array.
            slave_clstr=()

            # Iterate \$j from \$i + 1 to \$seg_no.
            for (( j=i+1; j<seg_no; j++ ))
                do
                    tclstr="\$indiv-\$j.vs.\$i.clstr"

                    # If the tclstr file exists and is not empty, add it to the slave_clstr array.
                    if [ -s "\$tclstr" ]; then
                        slave_clstr+=("\$tclstr")
                    else
                        echo "No file \$tclstr"
                       exit 1
                    fi
                done

        # If the slave_clstr array is not empty, merge the master and slave cluster files.
        if [ \${#slave_clstr[@]} -gt 0 ]; then
            tclstrs=\$(echo "\${slave_clstr[@]}" | tr ' ' ' ')
            echo "\$i\t\$master_clstr \$tclstrs" >> \$out_clstr
        else
            # Otherwise, directly append the contents of master_clstr to the output file.
            echo "\$i\t\$master_clstr" >> \$out_clstr
        fi
    done


    cat "\${reps[@]}" > "\$out"

    """

    stub:
    """
    touch clstr_order.txt
    touch unique.fa
    """

}
