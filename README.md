Project by Nate Gaylinn:
[github](https://github.com/ngaylinn),
[email](mailto:nate.gaylinn@gmail.com),
[blog](https://thinkingwithnate.wordpress.com/)

# Overview

This project is a demonstration of an experimental evolutionary algorithm
design (tentatively named an "epigenetic algorithm"). It's inspired by a simple
observation: a gene sequence alone does not describe an organism, because amino
acids don't have any inherent *meaning* in nature. Instead, the cell must
*interpret* the gene sequence in order to build and operate a body, and that
interpretation process was evolved along with the gene sequence. This project
mimics that structure, evolving programs to solve a `FitnessGoal` while also
evolving an interpreter for those programs. The result is an algorithm that can
search over a large space of possible phenotypes by learning which regions of
that search space are best for evolving solutions to each `FitnessGoal`.

This project is based on an earlier
[prototype](https://github.com/ngaylinn/epigenetic-gol-prototype). Since then,
the goals of the project have become clearer. The design is simpler, in that it
does fewer things, but more complicated, in that it does them in a more
optimized way.

To explore this project and its results in a more visual format, check out the
slide-deck summary
[here](https://docs.google.com/presentation/d/1ZIYj8Rg4xPHukQ-Glk9_IZTgimEYh5dtWsyXrbkX9Dc),
or the output from running this project
[here](https://github.com/ngaylinn/epigenetic-gol-v1-results).

This algorithm is not yet particularly useful. It's meant as an exploration of
a wild idea, to justify [further exploration](#future-work). It is now the
basis for a scientific paper (not yet published) exploring evolvability.

## Game of Life

This demo uses Conway's Game of Life (GOL) as a playground for evolutionary
algorithms. For those not familiar with the GOL, you may want to check out this
[introduction](https://conwaylife.com/wiki/Conway's_Game_of_Life) and try
playing with this [interactive demo](https://playgameoflife.com/).

The GOL was chosen for a few reasons. First and foremost, this project is meant
to evolve complex phenotypes. That is, it finds *programs* that generate
solutions, not just solutions themselves. This will be important for [future
work](#future-work), where this algorithm may be adapted to different sorts of
programming tasks. As computer programs go, a GOL simulation is about as simple
as they come. It takes no input, it's fully deterministic, has no dependencies,
and aside from the game world itself, it has no state or output. Yet, it's
still complex, often producing surprising emergent behaviors from those simple
ingredients.

Other reasons for using the GOL are that it's relatively well known and
produces nice visuals that make it clear how well the genetic algorithm is
performing. It's also easy to optimize with [parallel
execution](#inner-loop-design).

That said, doing cool things with the GOL is *not* a primary goal for this
project. If you're a GOL aficionado and would like to help make this work more
interesting and useful to your community, I'd love to hear from you. Please
[contact the author](mailto:nate.gaylinn@gmail.com) for possible
collaborations.

## Motivation

A major practical limitation of evolutionary algorithms is that they require an
expert programmer to optimize the algorithm design to fit their specific use
case. The programmer carefully decides what range of phenotypes to explore, and
an encoding scheme to describe them with a gene sequence. They often devise
clever mutation and crossover strategies to make the search process more
efficient. With luck and a great deal of trial and error, this sometimes
produces useful results. The problem with this is that it's a lot of hard work
by the programmer, and the results are limited by their creativity and
foresight. Why do that by hand, when life can do the same thing via evolution?

In biology, the mechanisms of *epigenetics* serve to manage the interpretation
and variation of the gene sequence. This project explores one way of doing the
same thing in an artificial context. Rather than fine-tuning the performance of
a single genetic algorithm designed for a single `FitnessGoal`, this algorithm
is more like searching over many possible genetic algorithms to find which one
gets the best traction on a particular `FitnessGoal`.

If successful, this might one day lead to genetic algorithms that are more like
deep learning foundation models. A powerful generic model could be built by
specialists, then automatically fine-tuned to a specific narrow task on demand.
It might also produce results that are more open-ended, flexible, and life-like
than traditional evolutionary algorithms.

In biology, the notion that life "evolves for evolvability" is still
controversial. The hope is this line of research might produce evidence in
support of that hypothesis. Although not very [biologically
realistic](#biological-realism), this project suggests that "evolving for
evolvability" is an effective strategy that fits naturally with the sort of
nested "program and interpreter" design seen in nature (for more on how the
cell resembles an interpreter, check out my [blog
series](https://thinkingwithnate.wordpress.com/2023/08/02/a-cells-eye-view-of-evolution/)
on the topic).

# Technical Overview

This project is a special kind of meta-evolutionary algorithm. It's basically
two evolutionary algorithms, one nested inside the other. The ["outer
loop"](#outer-loop-design) of the project evolves species&mdash;subsets of
possible organisms, each with a range of possible forms. A species is
represented by a `PhenotypeProgram`, which is used to turn a `Genotype` into a
phenotype (in this case, a GOL simulation). The ["inner
loop"](#inner-loop-design) holds the set of species constant and evolves a
population of organisms (each represented by a `Genotype`) for each
`PhenotypeProgram`. It uses a process of
[development](#genetics-and-development) to combine a `PhenotypeProgram` with a
`Genotype` to make a GOL simulation from each organism, then scores those
simulation based on one of several different `FitnessGoal`s. The fitness of a
species in the outer loop is determined by how effective the inner loop is at
evolving fit organisms using that `PhenotypeProgram`.

The outer loop is implemented in Python, mostly because the code is simpler,
more readable, and doesn't have to be high performance. It's purpose is to
randomize and breed `PhenotypeProgram`s, run the inner loop, and analyze the
results. The inner loop is written in C++ and is optimized for execution on a
CUDA-enabled GPU. It's main purpose is to run many GOL simulations in parallel.
It also evaluates the fitness of those simulations, and handles randomizing
and breeding of organism populations. The primary interface between these two
halves of the project are the `Simulator` class (which manages populations of
organisms on the GPU) and the `PhenotypeProgram` data structure (which the
outer loop evolves and the inner loop uses for development).

## Nested Evolution

The primary challenge with evolving a program to interpret gene sequences is
*stability*. A traditional evolutionary algorithm requires a stable and
relativel smooth [fitness
landscape](https://en.wikipedia.org/wiki/Fitness_landscape#:~:text=In%20evolutionary%20biology%2C%20fitness%20landscapes,often%20referred%20to%20as%20fitness)
in order to "get traction" on a problem. If the signal is too noisy or
unreliable, then it can't tell which mutations are an improvement, and
selection will no longer lead it in fruitful directions. Allowing the meaning
of the gene sequence to change while evolution is happening completely breaks
this. Even a small mutation to a `PhenotypeProgram` completely reshapes the
fitness landscape, making hill climbing impossible. Or, from another
perspective, changing the `PhenotypeProgram` breaks any `Genotype`s evolved for
it. To work around that problem, this project breaks evolution into two phases.
The outer loop continuously evolves a popoulation of `PhenotypeProgram`s
(species). The inner loop evolves a population of `Genotype`s (organisms), but
*not* continously. In each generation of the outer loop, the inner loop starts
over with a randomized organism population. While this makes nested evolution
possible, it's not [biologically realistic](#biological-realism), and requires
[special intervention](#genetics-and-development) for `Genotype`-level
innovations to persist as species evolve.

![Diagram of control flow in this nested-evolution
architecture](assets/execution-flow.png)

Normally, to evolve GOL simulations with interesting qualities, a programmer
would use their intuition to hand-design an evolutionary algorithm, using what
they know about GOL and a specific `FitnessGoal` to make the search process
efficient enough to work. A major [motivation](#motivation) for this project is
to avoid doing that, to search more broadly with less human bias. Instead, this
algorithm is optimized to search over the domain of monochrome bitmap images,
with no awareness of the GOL or any specific `FitnessGoal`. This project provides
many `FitnessGoal`s, chosen so that different evolutionary strategies are
needed to solve them. This shows that a epigenetic algorithm can find effective
evolutionary strategies to a variety of problems, all on its own.

An epigenetic algorithm produces two outputs: the best `PhenotypeProgram`
(species) and `Genotype` (organism) for a given `FitnessGoal`. These are
combined to make the organism's phenotype, which is a GOL simulation. The
fitness of that simulation is determined by counting which `Cell`s are `ALIVE`
and `DEAD` at different `Frame`s, looking for patterns that suit the
`FitnessGoal`s (in other words, they have properties that the programmer
decided were "desirable"). The fitness of a species is derived from those
organism fitness scores, but in an indirect way. Species are evolved for
*evolvability*, not fitness. As such, species fitness is determined by the
performance of the inner loop&mdash;how well fitness of the whole population
improved over many generations. This is approximated using a metric called the
"weighted median integral" (see
[`evolution.compute_species_fitness()`](evolution.py#L39))

Since the inner and outer loops of this project have different notions of
fitness, it's best to think of them as two separate evoluationary searches. The
inner loop searches for interesting GOL simulations within the search space
defined by a particular `PhenotypeProgram`. The outer loop searches for
`PhenotypeProgram`s, which represent subsets of the larger search space where
interesting GOL simulations can be found. In other words, an epigenetic
algorithm doesn't search for fit organisms, it searches for evolutionary
algorithms that produce fit organisms.

To get a sense for how different species specialize in different kinds of
phenotype, and for the range of different phenotypes explored by this project,
check out some of the examples and slide notes in the [overview
presentation](https://docs.google.com/presentation/d/1ZIYj8Rg4xPHukQ-Glk9_IZTgimEYh5dtWsyXrbkX9Dc/edit#slide=id.g28c7ff4855e_0_37).

## Genetics and Development

(this section covers [`kernel/phenotype_program.h`](kernel/phenotype_program.h) and
[`kernel/development.cu`](kernel/development.cu))

The core challenge of this project is to evolve a grammar for the gene sequence
of an evolutionary algorithm. Like any evolutionary algorithm, this project
must evolve a `Genotype` that effectively solves a FitnessGoal. What's
different is that the *meaning* of that `Genotype` is not fixed, but is
determined by a `PhenotypeProgram`, which gets evolved separately. The
`Genotype` itself, then, is mostly just a sequence of bits. For practical
reasons, those bits are organized into "genes" of two types: `Scalar` genes are
unsigned 32-bit ints, and `Stamp`s are 8x8 arrays of phenotype `Cell` values.
These genes take on meaning when they get bound to `Argument`s in a
`PhenotypeProgram`. `ScalarArgument`s are used to configure the behavior of an
`Operation` (translate *this many* pixels to the right), while `StampArgument`s
are used to encode portions of the phenotype literally (draw *this pattern*
here).

The `PhenotypeProgram`, is like a tree of primitive functions (`Operation`s)
that does nothing but transform data from the organism's gene sequence and
environment. For those familiar, it's a bit like the
[NEAT](https://towardsdatascience.com/neat-an-awesome-approach-to-neuroevolution-3eca5cc7930f)
algorithm for neuroevolution. This project uses two types of `Operation`s.
`DrawOperation`s actually write to the phenotype, usually by copying data from
a `Stamp` gene. `TransformOperation`s like `TRANSLATE`, `CROP`, `MIRROR` and
`COPY` don't change the phenotype, but instead warp its coordinate system. They
can either apply globally, which allows the `PhenotypeProgram` to position,
transform, and repeat the pattern in a `Stamp` gene across the GOL world, or
they apply to a `Stamp` gene, which can constrain what kinds of patterns will
evolve and get drawn. Overall, a `PhenotypeProgram` consists of one or more
`DrawOperation`s that layer on top of each other using some `ComposeMode`, with
zero or more `global_transform`s and `stamp_transform`s for each.

One important aspect of how `Argument` binding works is that each `Genotype`
actually has several genes of each type. Each `Argument` gets bound to one gene
of the appropriate type, but multiple `Argument`s can bind to the same gene
value. This allows for useful interactions between `Argument`s and
`Operation`s. For instance, a `TRANSLATE` `Operation` might always shift the
phenotype *diagonally* if the row and column offset `Argument`s are bound to
the same value. Similarly, a `CROP` and a `COPY` `Operation` could combine to
take an N-`Cell`-wide strip of a `Stamp` gene and repeat it exactly N `Cell`s
to the left, ensuring the two `Operation`s stay in sync even if N evolves to
take on a different value.

As mentioned above, the [nested evolution](#nested-evolution) strategy used by
this project has the unfortunate effect of throwing away evolved `Genotype`s
each time a new generation of species is spawned. To work around that problem,
this project uses the notion of `gene_bias`, which is loosely inspired by the
natural process of [canalization](#biological-realism). This works by randomly
fixing the value of a gene to some common value from the previous generation of
organisms. While this prevents the gene value from evolving into something even
better in the inner loop, it also allows the outer loop to quickly constrain
the search space to variations of a pattern that has worked before. For now,
`gene_bias` always works in this simple same way, though there are many
possible variations to explore, such as having a soft preference for a gene
value rather than a fixed value, or intentionally injecting randomness.

## Inner-Loop Design

The nested evolution architecture for this project requires running a *lot* of
GOL simulations. To evolve a population of species requires running many
simulated generations of each one, and a species generation requires evolving a
population of organisms over several generations. Making things worse, there's
a lot of variability in an evolutionary algorithm, so this project runs each
process of ecolution multiple times and considers the median performance.
Overall, a single experiment in this project requires 1.4 billion simulated
lifetimes (50 species x 5 trials x 150 generations x 50 organisms x 5 trials x
150 generations).

To make that feasible, the inner loop of this project runs on a CUDA-enabled
NVIDIA GPU. The code has been optimized so that in a development environment
with a single RTX a5000 GPU it can compute ~300,000 lifetimes per second. This
means a single species evolution experiment takes about 78 minutes on average
(some `PhenotypeProgram`s take longer than others). Each experiment is broken
into 5 species trials that take only 15-20 minutes to complete. The trials can
be run independently and in any order, making it possible to get partial
results for an experiment pretty quickly. Running further trials just tightens
the error bars on the results.

In GPU programming, most performance problems come from context switching and
memory transfers. To minimize these costs, this project runs 5 trials of 2500
organisms for a single generation of 100 simulated time steps as one batch
operation. Further generations reuse the same memory allocations. This is
implemented using the
[`Simulator`](kernel/simulator.cu)
class, which manages all the GPU-side memory for this project and orchestrates
all the function calls needed to evolve organisms and retrieve the results for
analysis. It's relatively cheap to construct a new `Simulator` object, but that
should still be avoided in an evolution loop that runs thousands of times. For
this project, each species trial gets its own `Simulator`, which performs well
and avoids a class of bugs where one trial might accidentally influence
another.

To minimize memory transfers between the GPU and host computer, `Genotype` data
for the whole population stays GPU resident the whole time. That means
randomizing and propagating `Genotype`s is performed on the GPU, not the CPU.
This operation has much less parallelism than running GOL simulations, so it
doesn't fully utilize the GPU device, but it's still more efficienty than
shuttling `Genotype` data back and forth across the PCI bus on every organism
generation. For randomness, this part of the project uses the cuRAND library.
Breeding organisms is a relatively simple process.
Aligning `Genotype`s for crossover is trivial, since each one is just an array
of genes, and only organisms of the same species are ever bred together. When
two `Scalar` genes cross over, the value is taken from one of the parents
randomly. When two `Stamp` genes cross over, half of the `Stamp` pattern will
come from one parent, half from the other. (both randomization and breeding can
be found in [`kernel/reproduction.cu`](kernel/reproduction.cu))

Selection (see [`kernel/selection.cu`](kernel/selection.cu)) is implemented using an algorithm known as [Stochastic Universal
Sampling](https://en.wikipedia.org/wiki/Stochastic_universal_sampling). This
technique will randomly make N selections from a population, where each
individual is chosen in proportion to its fitness and may be chosen zero, one,
or many times. Unlike other selection algorithms, however, this one avoids
unlikely (but possible and unhelpful) outcomes, like selecting the same
individual every time. When used for organism breeding, this operation is
performed on the GPU, simultaneously computing selections for 250 populations
of 50 organisms in one go. To avoid implementing this algorithm twice, the
`select` method is also exposed as an API for outer loop to use on the CPU to
perform selection on a single population of 50 species. This means the code is
written once, but compiled for both the GPU and CPU devices.

Most of the inner loop's time is spent running GOL simulations (see
[`kernel/gol_simulation.cu`](kernel/gol_simulation.cu)).
This begins by building the phenotype, which is the GOL world for the first
time step of the organism's simulated lifetime. This involves running the
species-specific `PhenotypeProgram`s, passing in each organism's `Genotype`, as
described [above](#genetics-and-development). Next, the GOL simulation is run
for 100 time steps and fitness scores are computed. Since the GPU has very
limited high-performance memory, actually recording each step of the GOL
simulation slows things down dramatically and is only done on demand. Instead,
each `Frame` usually overwrites the memory used to compute the previous
`Frame`, which makes it impossible to get a view of the organism's full
lifetime when computing fitness. This project uses a `FitnessObserver` class to
work around this limitation (see
[`kernel/fitness.cu`](kernel/fitness.cu)).
The `FitnessObserver` considers each step of the simulated life time as it's
computed and records whatever data it needs to compute an overall fitness score
when the simulation ends.

The interface between the inner loop (written in C++) and the outer loop
(written in Python) is a Python module named `kernel` and implemented in
[`kernel/python_module.cc`](kernel/python_module.cc). It's implemented using
[pybind11](https://github.com/pybind/pybind11), which provides a simple macro
language for mapping identifiers across the two languages. The `kernel` module
most be [compiled](#execution) before it can be used in Python code. The
`kernel` module consists of the `Simulator` class, the `PhenotypeProgram` and
`Genotype` data structures (C++ structs that appear in Python as
[Numpy](https://numpy.org/) dtypes), the core constants configuring inner loop
behavior, and a few unbatched operations for use with small numbers of
organisms / species or for testing.

## Outer-Loop Design

The outer loop of this project is responsible for evolving `PhenotypeProgram`s.
It starts by randomly generating a population of 50 species (known as a
`Clade`) (see [`evolution.py`](evolution.py)). This is done by starting with a single individual that has a minimal
`PhenotypeProgram` (draw a single `Stamp` with no transformations), and
systematically generating a diverse set of mutant clones of that individual to
fill out the rest of the population. It then invokes the inner loop to evolve
organisms for all these species, and computes species fitness (see
[above](#nested-evolution)) from the results. Cross breeding
`PhenotypeProgram`s is more complicated than `Genotype`s because their
structure is more complex and open-ended. This project uses the concept of
"innovation numbers" (borrowed from the [NEAT
algorithm](https://nn.cs.utexas.edu/downloads/papers/stanley.ec02.pdf)) to
align tree structures and produce hybrid trees (breeding for
`PhenotypeProgram`s is implemented in
[`phenotype_program.py`](phenotype_program.py)).

Since the concept of an epigenetic algorithm is novel, and it's not yet
entirely clear what it's good for or how to implement well, this project
explores a few variations of the algorithm to see what's important. These
variations are managed using the [`Constraints`](phenotype_program.py#L187)
dataclass.
It consists of three independent boolean variables: whether to use gene bias /
simulated canalization or not, whether to allow multiple draw operations or
just one, and whether to allow transformations on `Stamp` genes or only global
transformations (see [above](#genetics-and-development)). Each experiment has a
fixed set of constraints and a fixed `FitnessGoal`, and this project runs many
experiments to see how the different constraints affect performance across many
`FitnessGoal`s. The code for actually running an experiment and collecting data
can be found in
[`experiments.py`](experiments.py).

The main entry points for this project are
[`run_experiments.py`](run_experiments.py)
and
[`visualize_results.py`](visualize_results.py).
The first of these takes a batch of experiments and runs them all, one trial at
a time. In order to produce partial results quickly, it uses a breadth-first
approach, running the first trial of all experiments before moving on to the
second trial, the third, and so on. When each trial completes,
`run_experiments.py` automatically invokes `visualize_results.py` in a separate
process to analyze the results and output a set of user-friendly summaries,
including data tables, simulation videos, and charts. The
`visualize_results.py` script can also be run independently, to iterate on the
analysis and visualization process without re-running the evolutionary
experiments.

The outer loop also provides a few development tools.
[`benchmark.py`](benchmark.py)
is used for tuning performance of the inner loop. It basically just invokes
that code a few times to get an accurate runtime measurement, then compares
that to prior runs.
[`trace_species_evolution`](trace_species_evolution.py)
is for debugging and fine-tuning the species evolution process. It runs the
first few generations of a single experiment, generating an in-depth summary of
the evolutionary algorithm results from the inner-loop. This helps visualize
what each `PhenotypeProgram` actually does, and how that influences the
next generation. There's also a suite of
[`tests`](https://github.com/ngaylinn/epigenetic-gol-v1-wip/tree/main/tests)
which mostly serve to document and validate the behavior of the C++ code of the
inner loop (it's tested in Python for convenience and personal preference). The
outer loop is relatively simple and easy to verify, so
[`test_phenotype_program.py`](tests/test_phenotype_program.py)
is the only test for the Python code.

## Execution

This project requires a CUDA-enabled GPU to run. It was developed on a machine
with a single RTX a5000 GPU ([Compute Capability
8.6](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capabilities).
It could likely be run on other GPU devices without trouble, but this may
require some mucking around with the thread and block constants in
`environment.h`, the convenience functions in `cuda_utils.h`, or the build
configuration in `CMakeLists.txt`. It also has pretty large memory
requirements. To get it running on a GPU with less memory, try using smaller
population sizes or fewer trials / generations. It's also possible to remove
the `Video`s in the `Simulator::DeviceAllocations` data structure, if the
`ENTROPY` `FitnessGoal` is removed.

This project is built using CMake. It requires Python 3.10 (other version
numbers may also work), [pybind11](https://github.com/pybind/pybind11), the
standard CUDA development libraries, and the nvcomp library. It also depends on
several common Python modules, which can be installed with:

```bash
# For managing and analyzing experiment data
pip install numpy
pip install pandas
pip install scipy

# For making charts and data visualizations
pip install matplotlib
pip install seaborn

# For importing and exporting animated gifs
pip install Pillow

# For rendering progress bars on the command line
pip install tqdm
```

Before any of the Python code will work, and after any changes to the C++ code
the kernel module must be compiled. To setup CMake to build the kernel module,
run `cmake kernel`. After that, the kernel module can be built by running
`make kernel`. You can also run `make cpp_demo` to build a standalone C++
binary (`cpp_demo`) that runs the inner loop through its paces without any
Python code. This is primarily for analyzing performance of the C++ code using
Nvidia's CUDA debug tools, where the Python code would add nothing and could
possibly obscure things.

This project offers several Python scripts:
* `python3 run_experiments.py`: Runs all experiments found in
  [`experiments.experiment_list`](experiments.py#L321) and invokes
  `visualize_results.py` to render the results. Takes no arguments. This script
  can take many hours or even days to complete, depending on the list of
  experiments. Progress is saved after every trial (~20 minutes on the original
  development machine). Typing Ctrl-C once will gracefully exit after the
  current trial completes, and typing it a second time does a forced quit,
  losing any unsaved progress.
* `python3 visualize_results.py`: Reads experiment data produced by
  `run_experiments.py` and uses it to generate data tables, simulation videos,
  and charts summarizing the results. By default, this script renders outputs
  for any experiments with new data, but when passed the `--rebuild` argument,
  it will re-render outputs for all computed experiment data. Running this
  script before any experiments have run has no effect.
* `python3 benchmark.py`: Runs the inner loop through its paces, measures how
  long it takes, and compares that runtime against prior runs. Output is
  noramally just a data table printed to the console, but pass the `--chart`
  argument if you'd like a visualization of the past few runs.
* `python3 trace_species_evolution.py`: Runs a few generations of species
  evolution and produces an in-depth summary of each. Takes no arguments. To
  configure species evolution, just modify the constants at the top of the
  `main` function.
* `python3 -m unittest discover`: Runs all unit tests and reports any failures.
  These tests are mostly just a sanity check that everything in the inner loop
  is working properly. They also generate a weath of images in directories
  named `tests/test_*` which serve to visualize the outputs of various kernel
  operations and document what they do.

# Biological Realism

This project is inspired by cellular biology, but is not intended to model or
simulate a cell in a realistic way. Its main purpose is to explore the idea of
using an evolved program to interpet an evolved gene sequence. This idea is
taken from real life, but implemented in a way that's very differnet from how
cells do it.

Perhaps the least realistic part of the project is the [nested
evolution](#nested-evolution) architecture. Obviously, nature does not evolve
species in discrete generations, spawning a new randomized population of
organisms each time. It would be interesting to explore other, more natural,
continuous ways of doing this, but solving that problem wasn't necessary to
study evolved gene sequence interpreters.

The division of `PhenotypeProgram`, `Genotype`, and `make_phenotype` is also
pretty unrealistic. In nature, these concepts are tangled together into one
complex system, consisting of DNA molecules and a whole host of associated
protein complexes and cellular processes. It's also pretty weird that
development happens all at once when an organism is spawned, and then the rest
of the organism's life is computed from that phenotype in a fully deterministic
way. Development and living life are more dynamic, continuous, and interactive
for a natural organism. Finally, the term "development" is most often applied
to multicellular organisms that assemble and grow a complex body, but this
project operates more like a single cell directing its own growth. It does not
involve collaboration of many independent cells, each consulting their own gene
sequence, with cell-type differentiation, and so on.

The design of the `PhenotypeProgram` and `Genotype` bear very little
resemblence to their natural counterparts. Some aspects of the
`PhenotypeProgram` resemble things like methylation or hox genes, but only
obliquely.
[Methylation](https://en.wikipedia.org/wiki/DNA_methylation) is a natural
system used to annotate a gene sequence to change how it gets read and copied
by the cell. [Hox genes](https://en.wikipedia.org/wiki/Hox_gene) serve as a
general purpose body-building mechanism shared by all animals. The `Genotype`s
division into `Stamp` and `Scalar` genes is a bit like [protein-coding and
regulatory
genes](https://en.wikipedia.org/wiki/Human_genome#Coding_vs._noncoding_DNA) in
a biological cell. Resemblance to these things was *not* a goal of this
project, but emerged naturally from the requirements of the data structure,
which is interesting. It's also significant that `PhenotypeProgram`s do not
invoke genes in sequential order, but instead use random-access lookup. This
was intentionally borrowed from nature to make the results of development more
life-like, but the way its implemented is quite different from nature's [gene
regulatory networks](https://en.wikipedia.org/wiki/Gene_regulatory_network).

The concept of "gene bias" (see [above](#genetics-and-development)) was
inspired by
[canalization](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5217179/) in
nature. This is the phenomenon where an organism's most critical behaviors
become supported by many intertwined genetic pathways, making it so a failure
in any one of those pathways won't disrupt proper functioning. Conceptually,
canalization constrains what phenotypes are possible, making certain genetic
variations equivalent so mutations in those places have no effect. Gene bias
works very differently, but serves the same purpose of constraining the search
space to just phenotypes that include certain critical behaviors established by
prior generations.

Although it's true of pretty much *every* evolutionary algorithm, it's worth
noting that the way fitness and selection work in this project is profoundly
undrealistic. Nature doesn't have any sort of explicit fitness function or
authoritative judge of fitness. Organisms merely survive (or not) and reproduce
(or not) selecting a mate (or not) as they please.

Similarly, organisms normally compete and collaborate with each other in a
shared physical space, but in this project each individual is isolated. Their
environment consists of nothing but their phenotype, and the only interaction
they have with each other is at reproduction tie and micro-managed by the
evolutionary framework.
