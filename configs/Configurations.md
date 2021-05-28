Configurations for Ensembles
============================

Here are configuration files which define ensembles that shall be started during a simulation.
The configuration consists of a list of "segments" where each segment indicates an ensemble start.
The variable "when_method" can have the following values
- repeated_time: Create a new ensemble every x seconds, where x is defined via "when_value"
- impact_change: Create a new ensemble when the parameter with the highest impact to the model state variable defined via "when_sens" changes to a different parameter
- sign_flip: Create a new ensemble when the sign of the sensitivity of the model state variable defined via "when_sens" to the parameter defined via "when_name" is changed
- value: Create a new ensemble when the value given via "when_value" for a variable given via "when_name" is reached
The parameter "amount" defines the number of ensemble members. "when_counter" can be used to create multiple ensembles whenever the criterion specified via "when_method" is satisfied. This does not apply for "repeated_time". \
Each segment has a list "params", where the perturbed parameters for the ensemble are defined. Each param has the following entries:
- name: Name of the parameter to perturb
- sigma: The width of the possiblity distribution for the perturbance
- simga_perc: The width of the possiblity distribution for the perturbance as percentage of the parameter's absolute value
- type: Some parameters are defined for different hydrometeor types, which can be specified here.
- rand_func: Choose the distibution for the perturbance. Currently, only "uniform" is supported
Each ensemble can spawn new ensembles for all methods except "repeated_time" with the current implementation.