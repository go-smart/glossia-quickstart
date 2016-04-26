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

h3. Running without Glossia

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
appear in ``output/run/``. The code outputs most run-time numerics into the ``results/`` subdirectory.

h3. Running with local Glossia

Once you have <https://github.com/go-smart/glossia-server-side/> up and running, from the previous layout, run
    .. code-block::

        glot launch -i input/tumour.stl -i input/vessel1.stl -i input/vessel2.stl -i input/organ.stl --tmp-directory ${GLOSSIA_SERVER_SIDE_LOCATION}/transferrer original.xml code/*

This is actually a similar command to the launch above, but as ``glot`` does not process the ``original.xml`` definition,
(it leaves this to be done by Glossia) the user must specify the locations of any needed input surfaces. Moreover,
the tmp directory is given to allow exchange of files without a Glossia-side file server.

The results may be retrieved, as usual for Glossia, with the command,
    .. code-block::

        glot results ${SIMULATIONID}

If successful, this will output a TGZ bundle of simulation results in the local directory. For more diagnostic
data, check the `glot documentation <https://go-smart.github.io/glot>` for ``glot diagnostic``.
