
# Transparent Barrier Model for Spatial Gaussian Fields

## Abstract

Spatial Gaussian fields (SGFs) are essential tools in spatial and spatio-temporal modeling, with the Matérn model frequently applied due to its flexibility and computational advantages. However, standard SGFs assume stationarity and isotropy, conditions that become unrealistic in environments with physical barriers. The Barrier model addresses this limitation but assumes that barriers are fully impermeable, a significant restriction in real-world applications. To resolve this, we propose the Transparent Barrier model, an innovative approach incorporating barriers with varying permeability, thus providing greater modeling flexibility and realism. We illustrate the utility of this model with marine megafauna distribution, specifically Dugongs (*Dugong dugon*), in the Red Sea, demonstrating its potential to accurately reflect complex spatial structures involving partially permeable barriers while maintaining computational efficiency.

## Introduction

Spatial Gaussian fields (SGFs) are widely used in modeling spatial and spatio-temporal phenomena, particularly in applications where residual spatial structures arise due to unmeasured covariates, spatial aggregation, or spatial noise. Among these, the Matérn model is a prominent choice due to its flexibility and diverse applications. Advances such as the integrated nested Laplace approximation (INLA) and stochastic partial differential equation (SPDE) approach also support the use of the Matérn model by enabling efficient Bayesian inference and computational feasibility.

One weakness of SGFs like the Matérn model is their reliance on assumptions of stationarity and isotropy. However, these assumptions become unrealistic in the presence of physical barriers or irregular features such as coastlines or islands. The Barrier model was proposed to address this, replacing Euclidean distance with a structure that removes paths crossing physical barriers. However, it assumes full impermeability.

Many real-world scenarios involve barriers with varying permeability. The Transparent Barrier model extends the original Barrier model by incorporating partial permeability. For instance, islands may act as impermeable barriers, while sandbanks or tidal flats may be partially permeable. Transparency is introduced as a parameter controlling permeability, allowing flexible modeling of spatial dependency.

## Motivation

Marine megafauna such as Dugongs play vital ecological roles, particularly in seagrass ecosystem maintenance. Conservation efforts for such species are essential for sustaining marine biodiversity. Species distribution models (SDMs) are central tools in this domain, though their reliance on complete environmental datasets can be limiting.

In our case study along the northern Saudi Arabian Red Sea coast—characterized by island chains and irregular barriers—environmental data is limited, and bathymetry is used not as a covariate but to define physical structure. This scenario demands a model that can account for spatial random effects and variable barrier permeability.

The Dugong sightings are incidental, collected from tourism and citizen science. We employ a Poisson process model for these spatial point data. Our modeling approach captures underlying spatial structure and supports conservation goals in a data-limited, complex marine environment.

## Model

The Transparent Barrier model builds on the SPDE-INLA framework by assigning different Matérn fields to different regions, based on their permeability. The spatial domain Ω is divided into the normal area (Ωₙ) and barrier regions (Ω_b1, ..., Ω_bl), each with its own range parameter.

The SPDE for the field µ(s) is:
```
µ(s) - ∇·(r(s)^2 / 8 ∇µ(s)) = r(s) √(π/2) σµ W(s)
```
where r(s) varies by region. The resulting system is solved with the finite element method over a triangulated mesh. This yields a sparse precision matrix Q, allowing scalable inference.

## Correlation

The transparent barrier permeability is determined by the barrier range parameter, defined as a fraction of the spatial range in the normal area. Higher range fraction values increase barrier permeability, thereby strengthening spatial connections across barriers. Figure 1 illustrates how barriers transition from impermeable to completely permeable as the range fraction increases. Specifically, the first row corresponds to a fully impermeable scenario, which can be modeled using the existing Barrier model, whereas the last row corresponds to the stationary model scenario with no barrier effects. To capture correlation structures for all intermediate scenarios between these two extremes, the Transparent Barrier model becomes necessary.

When encountering a barrier, the distribution of possible locations for a moving point is influenced by the barrier's permeability. Under the Barrier model, the probability of finding this point on the opposite side of the barrier is nearly non-existent, as shown in the top-right plot of Figure 1. However, if an alternative path around the barrier exists, such as a canal or gap, some correlation may still occur, as seen in the top-left plot of Figure 1. In contrast, under a stationary Gaussian random field, the spatial correlation structure remains unaffected by any barriers, preserving consistent correlations regardless of the barrier's presence or the point's location within the study area (bottom row of Figure 1). 

The Transparent Barrier model aims to distort the spatial correlation that the moving point would exhibit under a stationary scenario, but without reaching the extreme behavior of the fully impermeable scenario. This distortion, controlled by adjusting the fraction of the range parameter used in the barrier area, depends on the specific characteristics and nature of barriers encountered in real-world applications.

**Figure 1:** Correlation plots illustrating spatial dependence for two barrier configurations. The first configuration (left side) represents a study area divided by a barrier containing a canal in the middle, connecting the upper and lower sections of the normal area. The second configuration (right side) depicts a study area completely divided by a thin barrier. Rows correspond to different barrier permeability levels, expressed by range fractions: 0.01 (first row), 0.2 (second row), 0.5 (third row), 0.7 (fourth row), 0.8 (fifth row), and 1 (sixth row). Columns correspond to distinct reference points from which spatial correlation is measured across the study areas.

## Simulation

**Simulation Procedure and Accounting for Barriers and Non-stationarity**  
The simulation involved defining distinct SPDE models for normal and barrier regions, assigning each region a specific spatial range parameter. Barriers were modeled by reducing the spatial range parameter to values of 0.01, 0.2, 0.3, 0.4, 0.5, 0.7, and 0.8, effectively weakening spatial correlation across these features. The precision matrix derived from the mesh and SPDE models was used to simulate spatial fields from a Gaussian distribution. Observational noise was then added to emulate real-world scenarios, resulting in simulated spatial data that accurately reflect the underlying non-stationary processes due to physical barriers.

**Comparison of Stationary, Barrier, and Transparent Barrier Models**  
To evaluate how each model captures the underlying spatial structure, we compared the stationary model, the original Barrier model, and the proposed Transparent Barrier model using posterior summaries. This comparison focused on the posterior spatial field, the posterior marginal standard deviation, and the estimated posterior range.

Figure 2 shows simulation results for two spatial configurations. On the left, configuration 1 features a barrier with a canal connecting the upper and lower regions of the normal area. On the right, configuration 2 introduces a thin continuous barrier that fully separates the domain. Each setting tests the ability of the models to reflect spatial structure in the presence of varying barrier permeability.

The posterior spatial field plots reveal notable differences in how each model interprets spatial continuity across the barrier. The stationary model, which assumes homogeneity across space, fails to account for the blocking effect of the barriers—this is evident in the overly smooth fields that do not attenuate across the barrier regions. In contrast, the Barrier model captures the discontinuity well in configuration 2 but overly restricts connectivity in configuration 1, where some permeability through the canal is expected. The Transparent Barrier model, shown in the larger plots of Figure 2, offers a more balanced representation. It maintains spatial continuity through the canal in configuration 1 while still reducing correlation across the barrier in configuration 2, where no alternative path is available.

The posterior marginal standard deviation highlights uncertainty associated with spatial prediction. In the Barrier model, uncertainty increases near the barrier, especially in configuration 2, reflecting the model’s inability to borrow strength across regions. The stationary model exhibits more uniform uncertainty but fails to reflect the true separation introduced by the barrier. The Transparent Barrier model shows a gradient in uncertainty that mirrors the structure of the domain: lower uncertainty where permeability allows information flow (e.g., near the canal in configuration 1) and higher uncertainty across impermeable regions (e.g., central regions of configuration 2).

Finally, the posterior range estimates are shown in the density plots alongside the spatial fields. The Barrier model often underestimates the effective range, since correlations are severely truncated at the barrier. Conversely, the stationary model overestimates range by ignoring structural boundaries. The Transparent Barrier model achieves intermediate and spatially adaptive estimates of the range, flexibly adjusting to local permeability.

**Figure 2:** Simulation results for two spatial configurations. On the left, configuration 1 features a barrier with a central canal, and on the right, configuration 2 shows a thin continuous barrier dividing the domain. Large colored plots display the posterior spatial field under the Transparent Barrier model. Smaller plots represent the Barrier model and stationary model, enabling visual comparison. Distribution plots show the posterior range estimates for each model and configuration.
