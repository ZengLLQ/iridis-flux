#!/usr/bin/env julia

using IridisFlux

exit_code = IridisFlux.run_cli(ARGS)
exit(exit_code)
