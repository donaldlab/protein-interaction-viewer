## Protein Interaction Viewer

### Installation

Protein Interaction Viewer (PIV) is a PyMOL plugin and can be installed like any [standard PyMOL plugin][3].   **Important:** In order to use this plugin the programs [Probe][PROBE] and [Reduce][REDUCE] must be in the user's path and be named probe and reduce respectively.  PIV also has a dependency on the [`future`][FUTURE] python package.  The [`future`][FUTURE] package should be installed either globally or for the user profile.

1. From the PyMOL toolbar select `Plugin-->Manage Plugins-->Install`. Or for newer versions of PyMOL go to `Plugin-->Plugin Manager` and then select the `Install New Plugin` tab.
2. Choose the PIV plugin file `ProteinInteractionViewer.py`
3. The plugin will show up in PyMOL menu `Plugin-->ProteinInteractionViewer`

#### Windows

- PIV does not work on the older Windows versions of PyMOL. This seems to be a problem with the version of python that these versions use. However, PIV does work with the newer open source versions of PyMOL. Instructions to find open-source pre-compiled binaries and installation instructions can be found [here][1].

#### Mac OSX

- The standard Mac OS PyMOL application, MacPyMOL.app does not support plugins by default. If you rename the PyMOL application to "PyMOLX11Hybrid.app" per the instructions [here][2]. PyMOL plugins will be enabled.
- The python version that MacPyMOL uses does not appear to add the user's path set in `~/.bash_profile` to its path. So please make sure probe and reduce are present in the system path (i.e. `/usr/bin`)

### Usage

#### Add and Remove Hydrogens

To add or remove hydrogens select the `Edit H` tab in PIV. Select the PyMOL object under the `Selection` box and use the `Clear H` or `Add H` buttons to remove or add hydrogens respectively. If no new object name is given, the new object will have the original object's appended by an "_H". If you would like to replace the current object, then check the `Replace` checkbox.

#### Load Contact Dots

To load contact dots select the `Load Dots` tab in PIV. Select the two PyMOL selections/objects you want to calculate dots between in the `Sel1` and `Sel2` combo boxes. Name the dots object, and put any additional parameters you want to send to the probe program. Click the `Load Dots` button to load the dots. In PyMOL the dots are grouped so that you can individually turn on and off different contacts. If you want to calculate dots between one object select that object in the `Sel1` box and then mark the `Self` checkbox.

#### Rotate Side Chains

Select the `Side-chain Rotator` tab to rotate side-chains. To select a sidechain to rotate, *ctrl+right clk* on any atom in the amino acid you would like to rotate (Note PIV doesn't update until you reselect the plugin window). Listed rotamers can be selected to change the side chain to that conformation. The scrollbars on the right allow for custom changes to each chi angle (Note after adjusting a chi angle with the scroll bar, smaller changes can be done using the left and right arrow keys).

#### General Comments

If the contact dots are hidden in PyMOL (i.e. the `hide everything` command was used), they can be shown again with the command show `cgo, <dots>` where *dots* is the name of the contact dots group.


### Citation

See [LICENSE](LICENSE)

### Other
See [the description on our website][WEBSITE] for more information about what this software does.

[1]: http://www.pymolwiki.org/index.php/Windows_Install#Pre-compiled_PyMOL
[2]: http://www.pymolwiki.org/index.php/Plugins
[3]: http://www.pymolwiki.org/index.php/Plugins_Tutorial#Installing_Plugins
[REDUCE]: http://kinemage.biochem.duke.edu/software/reduce.php
[PROBE]: http://kinemage.biochem.duke.edu/software/probe.php
[WEBSITE]: https://www2.cs.duke.edu/donaldlab/software/proteinInteractionViewer/index.php
[FUTURE]: https://pypi.org/project/future/ 
