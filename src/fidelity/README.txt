The files in this directory implement a different model from the temporal model.
Eventually, this may be given its own repository, but for now, they are just under 
the temporal model.

For background, see temporal/notes/6_12_17prosel_math_model.tex and the corresponding pdf file.

6/16/17

Implemented the framework to run the propsel retention (fidelity) simulation model
over multiple parameters, and to recored the result in a csv file.  This is similar
to previous code with a Fidelity module, a fidelity_type.

Sample run:

[wright@pardosa fidelity]$ julia -p 4 -L Fidelity.jl run.jl examples/example2

Code is currently in  evotech/temporal/src/fidelity, but if this idea is pursued,
it will be moved to evotech/fidelity, and a separate github repository created.

See also temporal/notes/6_12_17prosel_math_model.pdf and the corresponding LaTex file.
