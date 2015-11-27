# compass-fit example #
  1. In Materials Studio, create a very simple molecule, and manually set the force field types to the ones that you wish to fit. (example - fitting a three-atom group c3"-n3mh-c3a named **trio.xsd**.) The charges of each atom should remain zero.
> > ![http://pymsi2lmp.googlecode.com/svn/wiki/cnc.jpg](http://pymsi2lmp.googlecode.com/svn/wiki/cnc.jpg)
  1. In Discover->Setup->Automation Tab, set
    * Calculate forcefield types = No
    * Calculate partial charges = No
    * Calculate charge groups = No
  1. Now export the molecule as InsightII molecule files (in the example this will create **trio.car** and **trio.mdf**) to the compass-fit working directory.  In this example, you will find it in **"pymsi2lmp/examples/trio"**.
  1. Now open Discover->Dynamics, set
    * Temperature = 1000K
    * Number of steps = 100
    * Frame output every = 5
    * Trajectory: save = Coordinates
    * <font color='red'>(Number of steps / Frame output every) must be greater than the number of unknown force field parameters.</font>
  1. Run the simulation and copy **"Project Directory\Documents\trio disco dynamics\trio.arc"** to the compass-fit working directory.
  1. Enable the animation toolbar in Materials Studio and double click on  the **trio.xtd** in **trio Disco Dynamics** to open it.  Hit the back button on the animation toolbar to select the final frame of the structure.  At the status bar on the bottom you should see something like Frame 20/20.
  1. With this last frame still open, open Discover->Setup->Energy tab, and click calculate.  This will create a file **"Project Directory\Documents\trio Disco Dynamics\trio Disco Energy\trio.out"**.  Copy this file to the compass-fit working directory.
  1. Now you should have four files in the pymsi2lmp working directory (**trio.mdf**, **trio.car**, **trio.arc**, and **trio.out**).
  * <font color='red'>Note: unless you have LAMMPS compiled on Windows, you will probably need to copy this folder from your Windows workstation to a Linux machine.  At this point we are done with Materials Studio.</font>