Create input file
=================

Overview
--------------
We use ``.yml`` format for input files which has a dictionary-like structure. It has the following main sections:
    - ``mode``: str (optional), one of ``simulation``, ``postprocess`` or (default) ``simulation_postprocess``
            Type of calculation to perform.

    - ``fields``: dict of dicts
            Parameters of background fields.

    - ``grid``: dict
            Grid parameters.

    - ``integrator``: dict
            Type of amplitude to calculate (total vacuum emission or linearized in the probe field).

    - ``performance``: dict
            Performance-related parameters.

    - ``postprocess``: dict
            Postprocessing parameters (which observables to calculate from the complex amplitudes).


Fields
--------------
This section of ``.yml`` file is constructed as ``{'field_1': {...}, 'field_2': {...}, ...}`` where each dictionary of field parameters has the ``field_type`` parameter and 
other parameters specific to the chosen field type. 

``field_type`` is constructed via the combination of the field name (``dipole``, ``paraxial_gaussian``, ``laguerre_gaussian``) 
and how it is simulated (``analytic``, ``maxwell``). For instance, to use the analytic expression of paraxial Gaussian for all time steps, choose ``paraxial_gaussian_analytic``;
to use the dipole pulse expression for initialization and use linear Maxwell equations for later time steps, choose ``dipole_maxwell``.

Currently there are two special fields that are supported only as an analytic expression: ``eb_inhomogeneity``, ``plane_wave``.

For the full list of available keywords, refer to ``quvac.field.ANALYTIC_FIELDS`` and ``quvac.field.SPATIAL_MODEL_FIELDS``.

+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|                   |                                                                         Field types                                                                          |
| Field parameters  +----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|                   |                      dipole                        |                   paraxial_gaussian                |                  laguerre_gaussian                 |
+===================+====================================================+====================================================+====================================================+
|   ``focus_x``     |                                                                     Focus location in space                                                                  |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|   ``focus_t``     |                                                                      Focus location in time                                                                  |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|      ``W``        |                                                                            Energy                                                                            |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|     ``lam``       |                                                                          Wavelength                                                                          |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|     ``tau``       |                                                                           Duration                                                                           |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|``theta``, ``phi`` |        Virtual dipole moment orientation           |                                          Optical axis orientation                                       |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|    ``beta``       |                       --                           |                                             Polarization angle                                          |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|   ``phase0``      |                       --                           |                                          Phase delay at the focus                                       |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|     ``w0``        |                       --                           |                                               Waist size                                                |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|    ``order``      |                       --                           |                                          Paraxial expansion order                                       |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|     ``E0``        |                       --                           |                                   Field amplitude (if W is not specified)                               |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|     ``p``         |                       --                           |                         --                         |      Radial index of the Laguerre-Gaussian mode    |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|     ``l``         |                       --                           |                         --                         |    Azimuthal index of the Laguerre-Gaussian mode   |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|  ``dipole_type``  |           ``electric`` or ``magnetic``             |                                                   --                                                    |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+
|   ``envelope``    |  Temporal envelope type (``plane`` or ``gauss``)   |                                                   --                                                    |
+-------------------+----------------------------------------------------+----------------------------------------------------+----------------------------------------------------+

Grid
--------------
Required keys are:
    -  ``mode`` : str
        Mode of grid creation (``dynamic`` or ``static``).

Keys for ``static`` mode:
    - ``box_xyz`` : tuple of float
        Box size for the spatial grid.

    - ``Nxyz`` : tuple of int
        Number of grid points along each spatial dimension.
    
    - ``Nt`` : int
        Number of temporal points.
    
    - ``box_t`` : float or tuple of float
        Time duration or start and end times for the temporal grid.

Keys for ``dynamic`` mode:
    - ``collision_geometry`` : str
        Specifies the collision geometry ('x', 'y', 'z').

    - ``transverse_factor`` : float
        Factor to scale the transverse size.
    
    - ``longitudinal_factor`` : float
        Factor to scale the longitudinal size.
    
    - ``time_factor`` : float
        Factor to scale the time duration.
    
    - ``spatial_resolution`` : float or list of float, optional
        Controls the spatial resolution.
    
    - ``time_resolution`` : float, optional
        Controls the temporal resolution.
    
    - ``ignore_idx`` : list of int, optional
        Indices of fields to ignore for dynamic grid creation.


Integrator (optional)
---------------------
Keys:
    - ``type``: str
        ``vacuum_emission`` (calculate the total vacuum emission amplitude) or ``vacuum_emission_channels`` (calculate the amplitude linearized in the probe field)
    - ``probe_pump_idx``: dict
        Indices of probe and pump fields.
            - ``probe``: list of int
                Indices of the probe field, by default [0].
            - ``pump``: list of int
                Indices of the pump field, by default [1].

Performance (optional)
----------------------
Keys:
    - ``precision``: str
        Numerical precision for calculations: ``float32`` or (by default) ``float64``.
    - ``nthreads``: int
        Number of threads to use for ``numexpr`` library, by default all available CPUs.
    - ``pyfftw_threads``: int
        Number of threads to use for ``pyfftw`` library, by default equal to ``nthreads``.
    - ``test_run``: bool
        Whether to do a test run to estimate the resources for the full calculation.
    - ``test_timesteps``: int
        Number of timesteps for a test run, by default 5.
    - ``use_wisdom``: bool
        Whether to use existing wisdom file for ``pyfftw`` planning.


Postprocessing
--------------
Relevant keys for the polarization-insensitive signals:
    - ``calculate_xyz_background`` : bool, optional
        Whether to calculate the background spectra on Cartesian grid, 
        by default False.
    - ``bgr_idx`` : int, optional
        Index of the background field, by default None.
    - ``calculate_spherical`` : bool, optional
        Whether to calculate the spectra on spherical grid,
        by default False.
    - ``spherical_params`` : dict, optional
        Parameters for the spherical grid, by default None.
    - ``calculate_discernible`` : bool, optional
        Whether to calculate the discernible signal, by default False.
    - ``discernibility`` : str, optional
        Type of discernibility, (by default) ``angular`` or ``spectral``.

Relevant keys for the polarization-sensitive signals:
    - ``perp_field_idx`` : int, optional
        Index of the perpendicular field, by default 1.
    - ``perp_type`` : str, optional
        Type of perpendicular polarization, ``optical_axis`` or ``local_axis``, by default None.
    - ``calculate_spherical`` : bool, optional
        Whether to calculate the spectra on spherical grid,
        by default False.
    - ``spherical_params`` : dict, optional
        Parameters for the spherical grid, by default None.
    - ``stokes`` : bool, optional
        Whether to calculate Stokes parameters, by default False.


