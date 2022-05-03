using Printf

bptype = "sbp"

base = "Dogon Longitudinal Study/BP_and_early_growth2"

function do_mediation()

    for v in ["HT", "WT"]
        for sex in ["female", "male"]
            for adj in [0, 1, 2, 3]
                # patterns 6 0 ==> -1->-1 versus 1->1
                c = `python mediation.py $v $bptype 1 $sex 6 0 $adj`
                run(c)
            end
        end
    end

end

function mediation()

    # Create a consolidated results file
    fn = "$bptype/mixed/dim_1/mediation.txt"
    tables = []
    out = open(fn, "w")
    for v in ["WT", "HT"]
        for sex in ["female", "male"]
	        for adj in [0, 1, 2, 3]
                f = @sprintf("%s/mixed/dim_1/%s_adj%d_%s_med.txt", bptype, v, adj, sex)
                s = readlines(f)

                # Copy the whole file to the target
                for line in s
                    write(out, line)
                    write(out, "\n")
                end

                # Save the tables in text form
                push!(tables, s[1])
                ii = findfirst(x -> startswith(x, "Childhood"), s)
                for i = 0:7
                    push!(tables, s[ii+i])
                end
                push!(tables, "")

                write(out, "\n\n\n")
            end
        end
    end
    close(out)

    # Save the mediation summary tables to a consolidated text file
    fn = replace(fn, "mediation.txt" => "mediation_tables.txt")
    out = open(fn, "w")
    for x in tables
        write(out, x)
        write(out, "\n")
    end
    close(out)
end

function upload()

    fn = "$bptype/mixed/dim_1/mediation.txt"
    c = `dropbox_uploader.sh upload $fn "$base/$bptype/mixed/dim_1"`
    println(c)
    run(c)

    fn = replace(fn, "mediation.txt" => "mediation_tables.txt")
    c = `dropbox_uploader.sh upload $fn "$base/$bptype/mixed/dim_1"`
    println(c)
    run(c)
end

if length(ARGS) != 1
    error("wrong argument")
end

if ARGS[1] == "mediation"
    do_mediation()
    mediation()
elseif ARGS[1] == "upload"
    upload()
else
    error("unrecognized command")
end
