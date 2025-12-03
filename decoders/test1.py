import neat

# XOR test cases: input -> expected output
xor_inputs = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]
xor_outputs = [(0.0,),    (1.0,),    (1.0,),    (0.0,)]

def eval_genomes(genomes, config):
    """Fitness function: evaluates how well each genome solves XOR."""
    for genome_id, genome in genomes:
        # Create a neural network from this genome
        net = neat.nn.FeedForwardNetwork.create(genome, config)

        # Start with perfect fitness, subtract error
        genome.fitness = 4.0

        # Test on all 4 XOR cases
        for xi, xo in zip(xor_inputs, xor_outputs):
            output = net.activate(xi)
            genome.fitness -= (output[0] - xo[0]) ** 2

# Load configuration
config = neat.Config(neat.DefaultGenome, neat.DefaultReproduction,
                     neat.DefaultSpeciesSet, neat.DefaultStagnation,
                     'config_memory')

# Create population
p = neat.Population(config)
p.add_reporter(neat.StdOutReporter(True))

# Run evolution for up to 300 generations
winner = p.run(eval_genomes, 300)

# Test the winner
print('\nBest genome:\n{!s}'.format(winner))