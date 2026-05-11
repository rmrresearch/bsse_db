# Load the structure
load /Users/ibrahims/Documents/Programming/undergrad_reasearch/ames_research_summer_2025/psi4/docstring_stuff/data/HF/optimal_tetramer.pdb

# Hide everything first
hide everything

# Show representations
show sticks
show spheres

# Set sizes
set sphere_scale, .2          
set stick_radius, 0.3 

# Apply color scheme
util.cbag                     

# Background and rendering settings
bg_color white
set antialias, 2
set ray_opaque_background, off

# Quality settings
set stick_quality, 16
set sphere_quality, 4
set ray_trace_mode, 1

# Position the view
zoom all
orient

# Optional: Save the session or image
# save HF10_visualization.pse
# ray 1200, 1200
png ~/Desktop/HF10_structure.png