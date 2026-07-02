# Clean the global display and isolate the interface pocket
hide everything, all
show cartoon, all
set cartoon_transparency, 0.5

# Color individual chains to optimize visual contrast
color cyan, chain A
color green, chain B
color magenta, chain C

# Select and render the target mutation hotspot site as a solid sphere
select target_site, chain C and resi 558
show spheres, target_site
color red, target_site

# Set publication standard rendering environment parameters
bg_color white
set ray_shadow, 0
set ambient, 0.4
