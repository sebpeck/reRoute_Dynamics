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

To mitigate anthropogenic climate change, the field of public transportation is steadily comitting to fully-electric vehicles. To effectively implement battery-electric buses, knowing the limits and dynamic energy needs of any given route is of particular value. Longitudinal Dynamics Models are a well established method for determining the energy requirements of a vehicle through straightforward physics - reRoute_Dynamics provides the tools to couple geospatial data and an LDM to generate models of route-specific drive cycles and their power requirements. This tool opens up the avenue for high-throughput insights into route optimization and battery aging. 

# Statement of Need

when testing battery-electric vehicles for their milage, the EPA standard [@epa_range_testing] is to run these EVs on a Dynamometer in accordance with one of several standard EPA driving cycles. The resulting milage is typically reduced by a factor of .7 to account for variability in driving conditions, driver aggresion, and HVAC use. However, this one-size-fits-all approach may not accurately encompass the variability in driving conditions, especially when it comes to larger vehicles being used for public transit. With frequent stops, sufficiently dynamic mass, and highly variable traffic and road conditions depending on geographic location, a standard drive cycle based solely on velocity has enough uncertainty it's worth taking a look at closer alternatives. 

![Standard Drive Cycle Velocity Profile.](Urban_Dynanometer.png)<br>
*Figure 1*: A standard EPA Drive Cycle [@epa_drive_cycle]

\autoref{fig:Realcycle} 

Longitudinal Dynamics Models are a well-established method of determining the energy use of a vehicle based on simple rigid-body physics and standard equations of motion. Previous work has shown that these models can effectively be used to model bus energy [@thesis],[@abdelaty] in software like simulink or python. However, these have been conventionally coupled with the aformentioned standard driving cycles. Instead, reRoute_Dynamics allows for the generation the drive cycle alongside the use of the LDM in accordance with provided geospatial data like elevation, path, stops, and signals. Additionally, while many previous papers handled acceleration as a simple value, modern busses are operated using acceleration profiles controlled by a governer. This set of tools allow for the implementation of such a profile in vehicle acceleration behavior, should the user desire that. Beyond this, the high customizability of parameters enables the option to generate Monte Carlo simulations of a given route to get a better understanding of the range of conditions the energy storage system may face. 

/autoref{fig:LDM}

While any one model created using this toolkit may not be able to encompass all of the intricate detail and variability of a bus's drive cycle, it can be used for rapid and reasonable estimation of a proposed route's requirements and estimated maintenance and charging downtime, and is modular enough that it can be used to suit all sorts of transit needs. For example, King County Metro in Washington is setting out to achieve net zero emissions by 2035. They have begun piloting fully electric busses, so a tool like this could be of particular use to them and other transit agencies looking to achieve the same goal. Much of this tool was developed with these agencies in mind, and many of the examples use data sampled from King County itself. 

/autoref{fig:Monte_Carlo}


# Features

At its core, reRoute_Dynamics has several submodules each used for handling different aspects of generating a drive cycle. Many of these have more in-depth tutorials contained within the repository, but to simplify it down:

1. Geography_Tools.py is used to handle geospatial data and format it such that it can be used by the other sub-modules.
This includes methods for getting the bearings between two points, the coordinates bounding a geospatial dataset, geographical point interpolation, reading raster elevation, among many others. The pinnacle of this is the Route class, which is used to contain the information for a given 'route' - its geometry, elevation, limits, stops, signals, and signs, as well as ways to save and load these in a dense .json format.

2. Object_Params.py contains the basic objects such as a bus, ESS, or trip that are used to contain parameters.
This is primarily used for modeling, saving, and loading object designs - like a bus' mass and frontal area, or a trip's expected ridership or wind bearing. This also includes an in-built method for generating acceleration profiles that are similar in nature to a route's standard. 

3. Physics_Engine.py handles the physics calculations of a route through an LDM, including braking, accelerating, and maintaining speed.
This is where many of the calculations around the Longitudinal Dynamics Model happen - how the wind and grade affect the drag on the vehicle, how velocity and motor power are determined by the external conditions.

4. Trip_Simulator.py is used to put all three of the above together and run a given Bus with an ESS for a trip along a defined route. 
This can act both as a standalone script, or a module - and it contains the driving logic that is used to determine what the next action the bus will take depending on the current conditions. When used as a standalone script, it uses a handful of saved files that are generated using the other modules, and outputs the results as a new file.

/autoref{fig:PackageInterplay}

# Limitations & Future Work

At the moment, this toolkit has been shown to make reasonable models of realistic bus driving cycles, and is capable of producing Monte Carlo simulations that can be used to provide insights into the variability in ESS load for a given route. While the accuracy and implications of these simulations are more pertinent to be adressed on their own, further development of the ESS calculations are warranted. This includes more dynamic modeling of cells, cell and module aging, and drivetrain mechanisms. Additionally, the geography handling and trip simulation take not-insignificant amounts of time, so the models produced are not reasonable for use in highly time-sensitive applications such as autonomous vehicles. 

# Citations


# Figures
![Standard Drive Cycle Velocity Profile.](Urban_Dynanometer.png)<br>
*Figure 1*: A standard EPA Drive Cycle
![Velocity Profile of King County Metro Route 49, collected by GPS.\label{fig:Realcycle}](r49cycle.png)

![Diagram of a Longitudinal Dynamic Model.\label{fig:LDM}]()
![Model comparison of velocity and acceleration distributions of a Monte Carlo simulation of KCM route 45.\label{fig:Monte_Carlo}]()
![Diagram of module interplay.\label{fig:PackageInterplay}]()
![\label{fig:}]()

# Acknowledgements

This project was originally based off of the "Route Dynamics" Package created by Dr. Erica Eggleton, so a special thanks to her for the original concepts and groundwork. 
Thank you to Dr. Daniel T. Schwartz, for advising me during the course of developing this, and to Dr. Dave Beck, whom helped advise me and laid the foundation for many of the skills used in the creation of this. 


