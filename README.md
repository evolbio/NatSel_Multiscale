# Circuits_01: Overview

Julia source code for manuscript:

Frank, S. A. 2024. Natural selection at multiple scales.

by Steven A. Frank, <https://stevefrank.org>.

The code and this file are on [GitHub](https://github.com/evolbio/NatSel_MultiScale.git).

This software provides the code to reproduce the figures in the manuscript.

# Getting started with the code

I outline the basic steps for installing Julia and running the code.

## Setup Julia

1. Install [Julia](https://julialang.org/) via the [juliaup](https://julialang.org/downloads/) program.
2. Download the source code for this project from this [GitHub](https://github.com/evolbio/NatSel_MultiScale.git) page. On that page, click on the **Code** button. The simplest option is to download and uncompress a zip archive. Alternatively, you will see options to clone the git repository if you want to do that.
3. I have tested the code only on MacOS but it should also run on Linux and Windows.

## Instantiate the environment for a project

1. In a terminal, go to the top directory of the code hierachy. In that directory is this file, README.md, LICENSE, and the subdirectory src/.

2. In this top directory for the code, run the project by typing in the terminal `julia --project=.`

4. Once you have started Julia by the prior instruction, you are in the [REPL interactive command line](https://docs.julialang.org/en/v1/stdlib/REPL/) that interacts with the Julia software. Alternatively, if you know what you are doing, you can use VS Code or other system.

5. Each project depends on a set of external Julia packages that must be installed. To install or work with packages, you need to switch the REPL to package mode. At the `julia>` prompt, before typing any other character, hit the `]` key. You should now see a prompt `(Project) pkg>`. If you started in the Reservoir project, the it will show `(MLS) pkg>`. 

6. Type `instantiate`, which will install all the required packages listed in the file Project.toml and all the related dependencies. That file lists the specific package versions that I used and tested, creating the same environment. If you change package versions by updating or otherwise, then you may encounter errors that prevent the code from running. You can also see which version of the Julia software I used by looking at the top of the Manifest.toml file in the project directory. Usually, any Julia version greater than or equal to that version will work. However, if you have a problem, you may want to try matching to the Julia version in the Manifest file. You can do that by using [juliaup](https://julialang.org/downloads/) to download and activate a specific version of Julia.

7. Once instantiate has finished, which may take several minutes, you should return to the main julia prompt. To do that, immediately after a new `(Project) pkg>` prompt, hit the Backspace key (or CTRL-C might work). You should then once again see the main`julia>` prompt. If you get stuck at any time, you can exit Julia or open a new terminal and start Julia again as above.

8. You should only instantiate once for each project. The next time you start the same project, it should be ready immediately without moving to the package prompt. If you move to another project for the first time, you do need to repeat this process. Each project has its own environment and has to be instantiated the first time the project is used.

# Running the code in a project

In the src/ subdirectory is a file, Run.jl. The code to reproduce one or more figures in the manuscript is in that Run.jl file. As described in Run.jl, copy and paste the relevant code to produce the figures.

# Understanding the code

There are not many comments in the code at present. To understand the code, look for the function call in Run.jl that calls a driver routine in the other source files.

# Notes

If you are going to modify the code, have a look at [Revise.jl](https://timholy.github.io/Revise.jl/stable/).

Occasionally there are specific directory location strings set within this source code that may not work on your computer. If you run into a problem, change the directory in the source code.
