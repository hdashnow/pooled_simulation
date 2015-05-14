

num_to_pool = [2, 4]

check_samples = {
    branch.num_samples = branch.name.toInteger()
    def all_inputs = "$inputs".split(" ").toList()
    forward(all_inputs[0..branch.num_samples-1])
}

print_samples = {
    exec "echo $input.bam"
}

run {
    num_to_pool * [ check_samples +  "%" * [print_samples] ]
}
