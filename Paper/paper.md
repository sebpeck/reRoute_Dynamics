---
title: 'reRoute_Dynamics: A Package for Battery-Electric-Bus Modeling'
tags:
  - Python
  - EVs
  - dynamics
  - Busses
  - Transit
authors:
  - name: Sebastian Sachet Peck
    orcid: 0000-0000-0000-0000 TODO
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Daniel T. Schwartz
    orcid: 0000-0000-0000-0000 TODO
    equal-contrib: False
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: David Beck
    orcid: 0000-0000-0000-0000 TODO
    equal-contrib: False
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Erica E. Eggleton
    orcid: 0000-0000-0000-0000 TODO
    equal-contrib: False
    affiliation: "1" # (Multiple affiliations must be quoted)

affiliations:
 - name: University of Washington, United States
   index: 1
   ror: 00hx57361 TODO
date: XX April 2025
bibliography: paper.bib
---

# Summary

To mitigate anthropogenic climate change, the field of public transportation is steadily comitting to fully-electric vehicles. To effectively implement battery-electric buses, knowing the limits and dynamic energy needs of any given route is of particular value to these transit agencies, as well as Electric Vehicle manufacturers and researchers. Standardized drive cycles and Longitudinal Dynamics Models are a well established methods for determining the energy requirements of a vehicle through straightforward physics - the reRoute_Dynamics package offers the the tools to couple geospatial data and an LDM to generate models of route-specific drive cycles and their power requirements. This opens up the avenue for high-throughput insights into route optimization and the dis-aggregation of use-specific battery aging modes.

# Statement of Need
When testing battery-electric vehicles for their milage, the EPA standard `[@epa_range_testing]` is to run these EVs on a Dynamometer in accordance with one of several standard EPA driving cycles. The resulting milage is typically reduced by a factor of .7 to account for variability in driving conditions, driver aggresion, and HVAC use. However, this one-size-fits-all approach may not accurately encompass the variability in driving conditions, especially when it comes to larger vehicles being used for public transit.

![Standard Drive Cycle Velocity Profile.](Urban_Dynanometer.png)<br>
*Figure 1*: EPA Standard for Urban Dynanometer Testing, adapted from:`[@epa_drive_cycle]`

![Velocity Profile of King County Metro Route 49, collected by GPS.](r49cycle.png)<br>
*Figure 2*: Velocity profile of King County Metro route 49, collected via GPS.

 With frequent stops, sufficiently dynamic mass, and highly variable traffic and road conditions depending on geographic location, a standard drive cycle based solely on velocity has enough uncertainty it's worth taking a closer look at alternative methods. 

![Longitudinal Dynamics Model.](LDM.png)<br>
*Figure 3*: Diagram of a Longitudinal Dynamics Model.

Longitudinal Dynamics Models are a well-established method of determining the energy use of a vehicle based on simple rigid-body physics and standard equations of motion. Previous work has shown that these models can effectively be used to model bus energy `{[@thesis],[@abdelaty]}` in software like simulink or python. However, these have been conventionally coupled with the aformentioned standard driving cycles. Instead, reRoute_Dynamics allows for the generation the drive cycle alongside the use of the LDM in accordance with provided geospatial data like elevation, path, stops, and signals. Additionally, because modern transit vehicles operate within margins set by a speed governer, `[@XDE35_Manual]` this set of tools allow for the implementation of such a profile in vehicle acceleration behavior. The high customizability of parameters enables the option to generate Monte Carlo simulations of a given route to get a better understanding of the range of conditions the energy storage system may face. 

![Monte Carlo Simulation output of Route 45.](r45monte.png)<br>
*Figure 4*: Selected values of a Monte Carlo Simulation performed on King County Metro Route 45 Outbound.

While any one model created using this toolkit may not be able to encompass all of the intricate detail and variability of a bus's drive cycle, it can be used for rapid and reasonable estimation of a proposed route's requirements. This information can extend into estimated maintenance and charging downtime, and is modular enough that it can be used to suit all sorts of transit needs. For example, King County Metro in Washington have begun piloting fully electric busses, so a tool like this could be of particular use to them and other transit agencies in how they approach their operations planning. 

![Example Load Profile of Route 45.](essload.png)<br>
*Figure 5*: Example ESS Load Profile of King County Metro Route 45. 

Likewise, for those developing electric vehicle technology and energy storage systems, it enables the ability to create realistic load cycles. These cycles can be used in the research of EV battery failure modes and help dis-aggregate the aging inherent to the Energy Storage System from the aging due to use conditions. 

# Features

At its core, reRoute_Dynamics has several modules, each used for handling different aspects of generating a drive cycle. Many of these have more in-depth tutorials contained within the repository, but to simplify it down:

### Geography_Tools.py 
Geography_Tools.py is used to handle geospatial data and format it such that it can be used by the other modules. This includes methods for getting the bearings between two points, the coordinates bounding a geospatial dataset, geographical point interpolation, and reading raster elevation. The pinnacle of this is the Route class, which is used to contain the information for a given 'route' - its geometry, elevation, speed limits, transit stops, signal lights, and signs, as well as ways to save and load these in a dense .json format. 

### Object_Params.py 
Object_Params.py contains the basic object classes such as a bus, ESS, or trip. These classes are used to contain parameters for use in the modeling. Values like a bus' mass and frontal area, or a trip's expected ridership or wind bearing are all contained here, and can be saved or loaded for re-use. When saved, these objects are stored in an easily editable text file, which can likewise be re-loaded using methods within Object_Params. As part of this, there is also a method for generating an 'acceleration profile'. While possible to model a bus according to a single set acceleration, many metro busses such as the XDE35 use governers and other systems to control the manner in which they accelerate `[@XDE35_Manual]`. So, the user is capable of providing their own custom profile, or generating one that behaves in a similar manner to the acceleration seen in the Braunschweig Drive Cycle `[@NREL_Drive_Cycle]`, here. 

### Physics_Engine.py 
Physics_Engine.py handles the primary physics calculations of a route through a Longitudinal Dynamic Model (LDM), including braking, accelerating, and maintaining speed. Critical to the model's response to external conditions is first how those conditions are calculated. As seen in Figure 3, a LDM is essentially a force balance. So, at any given point, the vehicle will be experiencing the external forces of wind ($F_w$) and road resistance ($F_r$), calculated as follows:
$$F_w = \frac{C_d *A_f*\rho}{2} *(v - v_w*\cos(\theta_{rw}))^2$$
$$F_r = (\sin(\theta_g) + \cos(\theta_g)*C_f) * g * m $$
Where $C_d$ and $C_f$ are the coefficients of vehcile drag and road friction, $A_f$ is frontal area, $\rho$ is air density, $v$ and $v_w$ are vehicle and wind speeds, and $\theta_{rw}$ and $\theta_{g}$ are relative wind angle and road grade angle (in radians). $g$ and $m$ are gravity and current vehicle mass. These are used to determine the acceleration the bus is experiencing at each point due to the environment, after dividing each by the total current mass of the vehicle.

The vehicle will also be experiencing 'internal acceleration', which is determined by wether the bus is braking, accelerating, or maintaining speed. When braking, the bus's decelelleration is calculated according to the following formula:

$$ a_{int} = (a_{br} * f_i * f_{br})$$

where $a_{br}$ is the maximum braking acceleration, $f_i$ is the inertial factor - the amount the bus's mass reduces the braking ability - and $f_{br}$ is the braking factor, how hard the driver is hitting the pedal. 

When maintaining speed, $a_{int}$ is set to match the current external acceleration of $\frac{F_w + F_r}{m}$, where possible. If this value exceeds what would otherwise be possible with the braking power of the bus, $P_{br-max} = -m*v*(a_{br}*f_i)$, or the motor power of the bus, the velocity is adjusted by using iterations of the standard kinematic equations. 

Acceleration is the most complicated of the three, given it is meant to adhere to a given velocity profile. By assuming the acceleration profile is operated via a speed governer, the profile is integrated and the closest velocity to the vehicle is used as the starting point. The profile is then stepped through, accelerating the vehicle by the requisite amount, so long as the motor power limit is not exceeded. If that power limit is hit, the bus's net internal acceleration will be zero. 

Each calculation performed will return a final velocity, a time change for that action, and a power use for that action. Depending on the operation, these will be calculated slightly differently - but generally, Power is calculated using $m*a*v$, where $a$ is the internal acceleration of the bus, prior to losses from the envoronment, $m$ is vehicle mass, and $v$ is vehicle velocity. 


### Trip_Simulator.py 
Trip_Simulator.py is used to put all three of the above together and run a given Bus with an ESS for a trip along a defined route. This can act both as a standalone script, or a module. Containing a single method, it combines a route that was generated using Geography_Tools.py, and a Bus, ESS, and Trip, and passes each point contained in the route through a system of logic bsaed on the conditions defined.  

![Flow Diagram of Process](image-2.png)
*Figure 6*: Flow diagram of the simulation process. 

This set of logic is used to determine which operation in the physics engine is used, which in turn informs the next point on the route.  

# Future Work

This toolkit has been shown to make reasonable models of realistic bus driving cycles, and is capable of producing Monte Carlo simulations that can provide insights into the variability in ESS load for a given route. While the accuracy and implications of these simulations are more pertinent to be adressed on their own, the toolkit is set up in such a way that its modularity allows for significantly more customization and fidelit - a model is only as good as the data provided to it. 

Future development should include better accounting for drivetrain mechanisms, more traffic and road condition effects, and improved ESS modeling. This could include things like an HVAC power use determined by speed and ambient temperature, or dynamic modeling of cells and module aging. Simulation of ridership changes and signal light timings are, at present, randomly generated according to values provided in the Trip object. To more accurately model individual routes in the future, these should be able to be provided by a user. Beyond this, the geography handling takes a not-insignificant amount of time, so the models produced are not reasonable for use in highly time-sensitive applications.

Nonetheless, this remains a powerful tool and scaffolding for further transit and energy research, and makes easier study and development of how different aspects of a route, a bus, and an ESS impact the effectiveness of a Battery Electric Vehicle. 


# Citations

-TO DO: See Paper.bib.

# Acknowledgements

This project was originally based off of the "Route Dynamics" Package created by Dr. Erica Eggleton, so a special thanks to her for the original concepts and groundwork. 
Thank you to Dr. Daniel T. Schwartz, for advising me during the course of developing this, and to Dr. Dave Beck, whom helped advise me and laid the foundation for many of the skills used in the creation of this. 


