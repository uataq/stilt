## Frequently Asked Questions

1. **Where do I go for help?**  
   - For general discussion and questions (not bug reports), please use the [Discussions forum](https://github.com/uataq/stilt/discussions)
   - If you want to file a bug report,
     - Ensure there are no related reports by [searching open issues](https://github.com/uataq/stilt/issues)
     - After verifying there are no related open issues, [open a new issue](https://github.com/uataq/stilt/issues/new)
2. **What resolution should I set for the footprint calculations?**
    - The resolution largely depends on the application.
    - For use with an emissions inventory, you generally want to match the resolution of the inventory.
    > The resolution should be in the units of the desired footprint projection (e.g., degrees for lat/lon, meters for UTM). The projection of the meteorological data is irrelevant.
3. **How do I verify the model ran correctly?**  
   - Check the simulation directory `<stilt_wd>/<output_wd>/by-id/<simulation_id>/` for the presence of output files:
     - `<simulation_id>_traj.rds`: Particle trajectory data
     - `<simulation_id>_foot.nc`: Footprint data

     as well as simulation log files:
        - `stilt.log`: STILT model log
        - `MESSAGE`: HYSPLIT messages
        - `WARNING`: HYSPLIT warnings
   - If `rslurm` was used to submit the job, check the Slurm output files in `<stilt_wd>/_rslurm*/` for any errors.
      > For `X-STILT`, check `XSTILT/rslurm_XSTILT`
   - To further diagnose issues, run the binary executable directly from the command line in the simulation directory: `./hycs_std`. This will print any errors to the terminal.