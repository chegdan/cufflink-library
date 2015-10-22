# How to Contribute #
This is an open source project that will improve with community support.  The project benefits from testing the code, using it to solve problems, helping with the development, and spreading the word to the OpenFOAM CFD community.  Here are some suggestions on how to contribute to the cufflink library.

## Testing ##
Test cases have been provided the distribution so that users can check that everything works fine and that simple performance tests can be run easily.  Some scripts to extract the times from the log files that are output from the solvers are also been provided.  If you have suggestions for improvement, please discuss on the [cufflink-users](http://groups.google.com/group/cufflink-users) group.  if you have bugs to report, please discuss them on the users group first before putting them on the issues page.

## Using the code ##
One of the most important ways to contribute back to the cufflink and OpenFOAM community is to use the code and share your experience of what works and what does not.  Use the code, figure out what problems are best for this library, and share.  We have found that problems that are heavy on the inner-iterations (linear system solvers) at each time-step are ideal for cufflink.  Using this approach of a `reTol` of zero will drive down continuity errors at each time-step.  If you use `relTol` to drive down your residuals and want more outer-iterations with a high residual of continuity errors, then cufflink may not be good for that problem.  Sharing these situations is key to algorithm development and reducing frustration in the long run.

## Developing ##
If you would like to contribute to the development of the code there are several ways:
  1. If you want to contribute new proconditioners to the library, then we can certainly add them and if they apply to [Cusp](http://code.google.com/p/cusp-library/), then it will definitely get passed on to that project.
  1. If you want to write a new linear system solver, that is not already included (or included in Cusp but not implemented here), please do so.  We could use some more asymmetric solvers.  We will be putting up a roadmap of what we are planning on doing in the near future.
  1. **Suggested improvements in algorithms** can be discussed in the [cufflink-users](http://groups.google.com/group/cufflink-users) group and will get included with the distribution.
  1. If you want to contribute ideas for the roadmap, discuss it on the
  1. For suggested features, you can go to the [Issues Log](http://code.google.com/p/cufflink-library/issues/list) or discuss it in the users group
  1. **If you contributed to cusp and want to add that feature to cufflink**, please let us know and we can let you add the code.