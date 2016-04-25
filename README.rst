Glossia FEniCS Example
======================

Requirements
------------

- docker
- Python3 (via pip3)
    - glot
    - glossia.comparator

Instructions
------------

- Clone this repository, containing toy surfaces and a basic FEniCS solver
- Run
    .. code-block::
    
        glot setup . code/*
        
  This prepares the directory with the necessary tools for running a Glossia container offline.
- Run
    .. code-block::
    
        sudo ./setup.sh
        sudo ./run.sh
    
This should download the necessary image from DockerHub and execute the simulation in the container. Output should
appear in `output/run/`.
