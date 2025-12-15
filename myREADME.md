
# Code description:

## 1. Functions
All the necessary functions to run the R code.

## 2. Steps
Toy example to illustrate the steps we would follow in a real life example to build the barrier model. Mainly getting the range fraction to plug in the barrier model.

## 3. two rb
Toy example when we have more than one barrier.

## 4. Examples
More geometry examples, with canals, etc.

Has:
+ Example 3
+ Example4.Rmd 
+ sq example.Rmd
+ sq^2 example.Rmd

## 5. Poisson

### Poisson data with toy mesh Example

Mesh here is built in the same way as previous examples (BTopic103). However, data is simulated as if it was a Poisson process (SPDE book).

### Poisson data with Poisson mesh Example.

The mesh built here is done according to what I would do in a real life application for point processes. The data is simulated as if it was a Poisson process (it's the same as previous example).

## 7. Discussion 

Has the discussions of the other files altogether 

## Untitled folder:
template witg JABES journal

article_clean.code: code directly to journal template
article_mesh: mesh with two barrier i.e. canal and mesh with just one. barrier in the middle. I'm trying to use the ´interior = seg´ from INLA

Untitled: is the actual article from the template
Template: is the actual template without my things

article_prev.code: code from Poisson example and maybe others?

Meeting 6 feb Havard Janet
I need too send them the intro and model, and meet with Janet next week
mesh and simulations: 
I need to check if I did the fraction or actually multiplied the prior by the prior*fraction.
Check whether the range is actually greater inside the barrier because im doing range^2 for the 0 transparency.
Check why this gives me the proper range estimation. I think if I do jjust the fraction the range is worse than then one obtained with the stationary, but maybe is because i'm doing the simulation with wrong range fraction too. Check
check parameter coverage not just posterior and confidence interval
expand domain cuz range rn is higuer than the domain if i squared it.
do i want to do kinda a continuous space so then the domain is actually 0 to 20 (or 10) like Haakon or do I actually want ti think about boundaries.
I think I should just extended to infinity
I don't need to do another mesh for point process then and only the simulation. Which would actually make my life a lot easier. This means doing it as Hakon so i have to extend the mesh and boundaries, use interior and no max.edge, etc. This part I can do for the applied case, i.e. Dugongs
Simulation to get the params back is just to prove our model i.e. code works, and then doing it with the stationary model just to prove is wrong. Which is a bit stupid cuz it's wrong by definition.
Do correlation plots... i'm not sure what they want me to do, but it's important 
I can actually change the fraction and square it or do whatever i want as long as i can justify it and it makes sense.
check whether i am estimating the range or the field. yo creo que si estimo el field but to do so i need the range so i'm actually doing both.

mesh should be in model description not simulation. I think everything should be on model description but why did it cchange with the seed if i was only changing the simulated data with seed right?
I need to check haakon's tutorials


Get data! from Raul, check sharks, Meghan and Cassie


article_mesh.pp: i'm trying to simulate and model poisson data but i'm not including it in the results because i can't do the true simulated field with inla.qsample for poisson














