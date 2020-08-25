using Printf

bptype = "sbp"

base = "box:Dogon Longitudinal Study/BP and early growth"

function do_mediation()

	for v in ["HT", "WT", "HAZ", "WAZ", "BAZ", "BMI"]
        for sex in ["all", "female", "male"]
            # -1->-1 versus 1->2
    	    c = `python mediation.py $v $bptype 1 dadbp $sex 0 6`
    	    run(c)
    	end
    end

end

function mediation()
    fn = "$bptype/mixed/dim_1/mediation.txt"
    tables = []
	out = open(fn, "w")
    write(out, "```\n")
	for v in ["HT", "WT", "BMI"]
		for s in ["female", "male", "all"]
			f = "$bptype/mixed/dim_1/$(v)_dadbp_med.txt"
			if s != "all"
				f = replace(f, "med"=>"$(s)_med")
			end
			s = readlines(f)

			# Copy the whole file to the target
			for line in s
                write(out, line)
                write(out, "\n")
            end

            # Save the tables in text form
            push!(tables, s[1])
            ii = findfirst(x->startswith(x, "Childhood"), s)
            for i in 0:6
                push!(tables, s[ii+i])
            end
            push!(tables, "")

            write(out, "\n\n\n")
		end
	end
    write(out, "```\n")
	close(out)

    # Save the combined results to one pdf file
	gn = replace(fn, ".txt"=>".pdf")
	c = `pandoc -o $gn $fn`
	println(c)
	run(c)
	c = `rclone copy $gn "$base/$bptype/mixed/dim_1"`
	run(c)
	println(c)

    # Save the mediation summary tables to a consolidated text file
    fn = replace(fn, "mediation.txt"=>"mediation_tables.txt")
    out = open(fn, "w")
    for x in tables
        write(out, x)
        write(out, "\n")
    end
    close(out)
	c = `rclone copy $fn "$base/$bptype/mixed/dim_1"`
	run(c)
	println(c)
end

if length(ARGS) != 1
	error("wrong argument")
end

if ARGS[1] == "mediation"
	mediation()
elseif ARGS[1] == "do_mediation"
	do_mediation()
else
	error("unrecognized command")
end
