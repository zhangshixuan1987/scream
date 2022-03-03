# SCREAM data structures

Here's a list of the most important data structures that are used in SCREAM.

* The [Field](concepts/field.md) class stack.
* The [AbstractGrid](concepts/grid.md) interface.
* The [AbstractRemapper](concepts/remap.md) interface.
* The [FieldManager and GridsManager](concepts/fm.md) classes.
* The [AtmosphereProcess](concepts/atm_proc.md) interface.
* The [IO](concepts/io.md) classes.
* The [SurfaceCoupling](concepts/sc.md) class.
* The [AtmosphereDriver](concepts/ad.md) class.

Talk also about where concrete instances of the abstract interfaces are.
E.g., talk about HommeDynamics, P3Microphysics, etc, or the AtmosphereProcessGroup.
Also talk about the very important PhysicsDynamicsRemapper.

Finally, find the place to talk about how scream hooks up to E3SM. Part of it
is the SurfaceCoupling, but there's more (all the content of `mct_coupling` folder).
Maybe AD, SC, and mct hooks should be in their own page?

In the rest of this page, give a bird-eye view of how the above classes are interconnected,
like "each atm proc receives input fields and compute output fields, stored in the
field manager. each field is defined on a particular grid, and grid remapping is handled
by the remappers. IO classes allow to read inputd (such as IC files), and write output
(including restart files). The AD stores everything, and is in charge of init-ing everything, calling
all routines in the right order during normal time-step, and cleanup at the end."

Details on pre-post conditions, impl details, etc., for each class/data struct will be given
in their own pages.
