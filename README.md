# HVplanar


Optional keyword arguments that controls FEMM settings
  hide - hide the FEMM window on opening and run in the background
  frequency - (inductance only) set the frequency of the current
  precision - set the desired precision for the solver
  min_angle - set the minimum anlge for the mesh generator
  material - add the properties of a material to the simulation environment (see xx for applying it to the specific region).
  segment_angle - change the maximum angle of arc segments. Affects mesh around arc segments. 
  
TransformerGeometry 


FancyTrack
  layer:int from 0 to layers on primary or secondary side
  track = the number counting from the inner side starting at 0, alternatively, 'inner' or 'outer' can be specified 
  side_h - 'inner' or 'outer'. Sets the side of the track to round on the r-axis
  side_v - 'high' or 'low'. sets the side of the transformer over or under the insulation layer
  rounding is 'single' or 'both' where 'single'. single rounds the corner closest to the insulation layer
  elongation in mm of track width (r direction).
 
Guard
    
