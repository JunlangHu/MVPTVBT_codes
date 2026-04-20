# Ultimate root finder
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons

def D2(m1, sigma, s):
    """D2: 6th degree polynomial"""
    coeffs = [
        -1,  # m1^6
        2*s,  # m1^5
        -s**2 - 3*sigma**2/2 + 2,  # m1^4
        2*s*sigma**2 - 2*s,  # m1^3
        -s**2*sigma**2/2 - 3*sigma**4/4 + 3*sigma**2/2 - 1,  # m1^2
        s*sigma**4/2 - s*sigma**2/2,  # m1^1
        -sigma**6/8 + sigma**4/4  # constant
    ]
    
    # Compute using Horner's method
    result = 0
    for coeff in coeffs:
        result = result * m1 + coeff
    return result

def D3(m1, sigma, s):
    """D3: 8th degree polynomial"""
    coeffs = [
        -3*sigma**2/2,  # m1^8
        4*s*sigma**2,  # m1^7
        -7*s**2*sigma**2/2 - 3*sigma**4 + 3*sigma**2,  # m1^6
        s**3*sigma**2 + 6*s*sigma**4 - 5*s*sigma**2,  # m1^5
        -7*s**2*sigma**4/2 + 2*s**2*sigma**2 - 9*sigma**6/4 + 3*sigma**4 - 3*sigma**2/2,  # m1^4
        s**3*sigma**4/2 + 3*s*sigma**6 - 11*s*sigma**4/4 + s*sigma**2,  # m1^3
        -7*s**2*sigma**6/8 + s**2*sigma**4/4 - 3*sigma**8/4 + 3*sigma**6/4,  # m1^2
        s*sigma**8/2 - s*sigma**6/8 - s*sigma**4/4,  # m1^1
        -s**2*sigma**6/8 - 3*sigma**10/32 + sigma**8/4 - sigma**6/8  # constant
    ]
    
    result = 0
    for coeff in coeffs:
        result = result * m1 + coeff
    return result

def D4(m1, sigma, s):
    """D4: 10th degree polynomial"""
    coeffs = [
        -s**2*sigma**4/2 - 9*sigma**4/4,  # m1^10
        3*s**3*sigma**4/2 + 27*s*sigma**4/4,  # m1^9
        -3*s**4*sigma**4/2 - 5*s**2*sigma**6/4 - 21*s**2*sigma**4/4 - 57*sigma**6/8 + 27*sigma**4/4,  # m1^8
        s**5*sigma**4/2 + 3*s**3*sigma**6 - 3*s**3*sigma**4/4 + 35*s*sigma**6/2 - 27*s*sigma**4/2,  # m1^7
        -9*s**4*sigma**6/4 + 3*s**4*sigma**4/2 - 5*s**2*sigma**8/4 - 45*s**2*sigma**6/4 + 21*s**2*sigma**4/4 - 141*sigma**8/16 + 57*sigma**6/4 - 27*sigma**4/4,  # m1^6
        s**5*sigma**6/2 + 9*s**3*sigma**8/4 + 3*s**3*sigma**4/2 + 33*s*sigma**8/2 - 83*s*sigma**6/4 + 27*s*sigma**4/4,  # m1^5
        -9*s**4*sigma**8/8 + 7*s**4*sigma**6/8 - 5*s**2*sigma**10/8 - 61*s**2*sigma**8/8 + 11*s**2*sigma**6/2 + s**2*sigma**4/2 - 171*sigma**10/32 + 159*sigma**8/16 - 57*sigma**6/8 + 9*sigma**4/4,  # m1^4
        s**5*sigma**8/8 + 3*s**3*sigma**10/4 + 5*s**3*sigma**8/16 + s**3*sigma**6/4 + 27*s*sigma**10/4 - 17*s*sigma**8/2 + 13*s*sigma**6/4,  # m1^3
        -3*s**4*sigma**10/16 - 5*s**2*sigma**12/32 - 27*s**2*sigma**10/16 + 7*s**2*sigma**8/8 - s**2*sigma**6/8 - 51*sigma**12/32 + 21*sigma**10/8 - 9*sigma**8/8,  # m1^2
        3*s**3*sigma**12/32 + s**3*sigma**10/16 - s**3*sigma**8/8 + 65*s*sigma**12/64 - 5*s*sigma**10/8 - 5*s*sigma**8/16,  # m1^1
        -s**4*sigma**10/32 - s**2*sigma**14/64 - s**2*sigma**12/64 - 5*s**2*sigma**10/32 - 3*sigma**14/16 + 5*sigma**12/16 - sigma**10/8  # constant
    ]
    
    result = 0
    for coeff in coeffs:
        result = result * m1 + coeff
    return result

def find_roots_in_range(func, m1_range, y_vals, tol=1e-10):
    """Find roots of a function in the given range"""
    roots = []
    
    for i in range(len(m1_range)-1):
        y1, y2 = y_vals[i], y_vals[i+1]
        
        # Check for sign changes
        if y1 * y2 <= 0:
            x1, x2 = m1_range[i], m1_range[i+1]
            
            # Linear interpolation for root approximation
            if y1 != y2:
                root_approx = x1 - y1 * (x2 - x1) / (y2 - y1)
                # Avoid duplicates
                if not any(abs(root_approx - r) < 1e-6 for r in roots):
                    roots.append(root_approx)
    
    return sorted(roots)

def find_positive_intersection_regions(m1_range, y1, y2, y3, threshold=1e-10):
    """Find regions where all three functions are positive"""
    positive_mask = (y1 > threshold) & (y2 > threshold) & (y3 > threshold)
    
    regions = []
    in_region = False
    start_idx = 0
    
    for i in range(len(positive_mask)):
        if positive_mask[i] and not in_region:
            in_region = True
            start_idx = i
        elif not positive_mask[i] and in_region:
            in_region = False
            regions.append((m1_range[start_idx], m1_range[i-1]))
    
    if in_region:
        regions.append((m1_range[start_idx], m1_range[-1]))
    
    return regions

# Create main figure
fig, ax = plt.subplots(figsize=(16, 10))
plt.subplots_adjust(bottom=0.25, left=0.1, right=0.95, top=0.95)

# Initial parameter values
sigma_initial = 1.0
s_initial = 0.5
x_min_initial = -0.7
x_max_initial = 1.7
resolution = 2000  # fixed resolution

# Initial x range
x_min, x_max = x_min_initial, x_max_initial
m1_initial = np.linspace(x_min, x_max, resolution)

# Compute initial values
y_D2 = np.array([D2(x, sigma_initial, s_initial) for x in m1_initial])
y_D3 = np.array([D3(x, sigma_initial, s_initial) for x in m1_initial])
y_D4 = np.array([D4(x, sigma_initial, s_initial) for x in m1_initial])

# Find initial roots
roots_D2 = find_roots_in_range(D2, m1_initial, y_D2)
roots_D3 = find_roots_in_range(D3, m1_initial, y_D3)
roots_D4 = find_roots_in_range(D4, m1_initial, y_D4)

# Plot initial curves
line_D2, = ax.plot(m1_initial, y_D2, 'b-', linewidth=2.5, alpha=0.8, label='D2 (6th degree)', visible=True)
line_D3, = ax.plot(m1_initial, y_D3, 'r-', linewidth=2.0, alpha=0.8, label='D3 (8th degree)', visible=True)
line_D4, = ax.plot(m1_initial, y_D4, 'g-', linewidth=1.5, alpha=0.8, label='D4 (10th degree)', visible=True)

# Plot root points - use different colors and markers
scatter_D2 = ax.scatter(roots_D2, np.zeros(len(roots_D2)), 
                       color='blue', s=120, marker='o', edgecolors='white', 
                       linewidth=2, label='D2 roots', zorder=10, alpha=0.9, visible=True)
scatter_D3 = ax.scatter(roots_D3, np.zeros(len(roots_D3)), 
                       color='red', s=100, marker='s', edgecolors='white',
                       linewidth=2, label='D3 roots', zorder=10, alpha=0.9, visible=True)
scatter_D4 = ax.scatter(roots_D4, np.zeros(len(roots_D4)), 
                       color='green', s=80, marker='^', edgecolors='white',
                       linewidth=2, label='D4 roots', zorder=10, alpha=0.9, visible=True)

# Add annotations for each root - store annotations for different curves
root_annotations_D2 = []
root_annotations_D3 = []
root_annotations_D4 = []

for root in roots_D2:
    ann = ax.annotate(f'D2: {root:.3f}', xy=(root, 0), xytext=(root, 0.05*(ax.get_ylim()[1]-ax.get_ylim()[0])),
                     arrowprops=dict(arrowstyle='->', color='blue', alpha=0.7),
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8),
                     fontsize=8, color='darkblue', ha='center')
    root_annotations_D2.append(ann)

for root in roots_D3:
    ann = ax.annotate(f'D3: {root:.3f}', xy=(root, 0), xytext=(root, 0.1*(ax.get_ylim()[1]-ax.get_ylim()[0])),
                     arrowprops=dict(arrowstyle='->', color='red', alpha=0.7),
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.8),
                     fontsize=8, color='darkred', ha='center')
    root_annotations_D3.append(ann)

for root in roots_D4:
    ann = ax.annotate(f'D4: {root:.3f}', xy=(root, 0), xytext=(root, 0.15*(ax.get_ylim()[1]-ax.get_ylim()[0])),
                     arrowprops=dict(arrowstyle='->', color='green', alpha=0.7),
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.8),
                     fontsize=8, color='darkgreen', ha='center')
    root_annotations_D4.append(ann)

# Find and plot positive overlapping regions
positive_regions = find_positive_intersection_regions(m1_initial, y_D2, y_D3, y_D4)
positive_patches = []
for region in positive_regions:
    patch = ax.axvspan(region[0], region[1], alpha=0.3, color='gold', 
                      label='All D>0' if region == positive_regions[0] else "")
    positive_patches.append(patch)

# Set initial axis limits
ax.set_xlim(x_min, x_max)
all_y = np.concatenate([y_D2, y_D3, y_D4])
y_min, y_max = np.min(all_y), np.max(all_y)
y_margin = 0.1 * (y_max - y_min) if y_max != y_min else 0.5
ax.set_ylim(y_min - y_margin, y_max + y_margin)

# Add reference lines
ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1.2)
ax.axvline(x=0, color='black', linestyle='-', alpha=0.5, linewidth=1.2)

# Labels and title
ax.set_xlabel('m1', fontsize=14, fontweight='bold')
ax.set_ylabel('D values', fontsize=14, fontweight='bold')
ax.set_title(f'D2-D4 Polynomials: σ={sigma_initial:.2f}, s={s_initial:.2f}, x∈[{x_min:.1f}, {x_max:.1f}]', 
             fontsize=16, fontweight='bold')
ax.grid(True, alpha=0.3, linestyle='--')

# Place legend at lower left
ax.legend(loc='lower left', fontsize=10, framealpha=0.9, bbox_to_anchor=(0.0, 0.0))

# Create slider area - placed at lower right
slider_height = 0.03
slider_spacing = 0.04
slider_start = 0.05

# Slider parameters
slider_width = 0.4  # slider width
slider_left = 0.5   # slider left start position

# σ slider - lower right
ax_sigma = plt.axes([slider_left, slider_start + 2*slider_spacing, slider_width, slider_height])
sigma_slider = Slider(ax=ax_sigma, label='σ value', valmin=0.01, valmax=1.5, 
                     valinit=sigma_initial, valstep=0.001, color='blue', alpha=0.7)

# s slider - lower right
ax_s = plt.axes([slider_left, slider_start + slider_spacing, slider_width, slider_height])
s_slider = Slider(ax=ax_s, label='s value', valmin=0.0, valmax=1.0, 
                 valinit=s_initial, valstep=0.01, color='red', alpha=0.7)

# x min slider - lower right
ax_xmin = plt.axes([slider_left, slider_start, slider_width, slider_height])
xmin_slider = Slider(ax=ax_xmin, label='x min', valmin=-3.0, valmax=2.0, 
                    valinit=x_min_initial, valstep=0.05, color='green', alpha=0.7)

# x max slider - lower right
ax_xmax = plt.axes([slider_left, slider_start - slider_spacing, slider_width, slider_height])
xmax_slider = Slider(ax=ax_xmax, label='x max', valmin=-2.0, valmax=3.0, 
                    valinit=x_max_initial, valstep=0.05, color='purple', alpha=0.7)

# Add curve control checkboxes - placed at lower left
checkbox_ax = plt.axes([0.05, 0.02, 0.15, 0.08])
check = CheckButtons(
    checkbox_ax, 
    ['Show D2', 'Show D3', 'Show D4'],
    [True, True, True]
)

# Add checkbox title
checkbox_ax.set_title('Curve Controls', fontsize=10, fontweight='bold', pad=5)

# Info text box placed at lower right (above sliders)
info_text = ax.text(0.98, 0.98, '', transform=ax.transAxes, fontsize=9,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, 
                             edgecolor='gray', pad=5))

def update_info_text(sigma, s, x_min_val, x_max_val, roots_D2, roots_D3, roots_D4, positive_regions):
    """Update info text"""
    # Basic parameters
    param_info = f"σ = {sigma:.3f}, s = {s:.3f}, x∈[{x_min_val:.2f}, {x_max_val:.2f}]"
    
    # Number and positions of roots
    root_info = f"D2 roots ({len(roots_D2)}): {[f'{r:.3f}' for r in roots_D2[:3]]}"
    if len(roots_D2) > 3:
        root_info += f" ... (+{len(roots_D2)-3} more)"
    
    root_info += f"\nD3 roots ({len(roots_D3)}): {[f'{r:.3f}' for r in roots_D3[:3]]}"
    if len(roots_D3) > 3:
        root_info += f" ... (+{len(roots_D3)-3} more)"
    
    root_info += f"\nD4 roots ({len(roots_D4)}): {[f'{r:.3f}' for r in roots_D4[:3]]}"
    if len(roots_D4) > 3:
        root_info += f" ... (+{len(roots_D4)-3} more)"
    
    # Positive overlapping regions
    if positive_regions:
        region_info = f"All D>0 regions: {len(positive_regions)}"
        for i, region in enumerate(positive_regions[:2]):  # show at most 2 regions
            width = region[1] - region[0]
            region_info += f"\n Region {i+1}: [{region[0]:.3f}, {region[1]:.3f}] (width={width:.3f})"
        
        if len(positive_regions) > 2:
            region_info += f"\n ... and {len(positive_regions)-2} more regions"
    else:
        region_info = "No regions where all D > 0"
    
    info_text.set_text(f'{param_info}\n\n{root_info}\n\n{region_info}')

def update_root_annotations(roots_D2, roots_D3, roots_D4):
    """Update root annotation positions - keep original positioning logic"""
    # Remove all existing annotations
    for ann in root_annotations_D2:
        ann.remove()
    for ann in root_annotations_D3:
        ann.remove()
    for ann in root_annotations_D4:
        ann.remove()
    
    root_annotations_D2.clear()
    root_annotations_D3.clear()
    root_annotations_D4.clear()
    
    # Get current y-axis range
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    
    # Add new annotations - keep original positioning logic
    offset_counter = 0
    
    for root in roots_D2:
        # Assign different vertical positions for each root
        y_offset = -0.05 * 2 * y_range + 0.02 * 2 * y_range * (offset_counter % 3)
        ann = ax.annotate(f'D2: {root:.3f}', xy=(root, 0), 
                         xytext=(root, y_offset),
                         arrowprops=dict(arrowstyle='->', color='blue', alpha=0.7, lw=1),
                         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.9, edgecolor='blue'),
                         fontsize=9, color='darkblue', ha='center', fontweight='bold')
        root_annotations_D2.append(ann)
        offset_counter += 1
    
    for root in roots_D3:
        y_offset = -0.10 * 2 * y_range + 0.02 * 2 * y_range * (offset_counter % 3)
        ann = ax.annotate(f'D3: {root:.3f}', xy=(root, 0), 
                         xytext=(root, y_offset),
                         arrowprops=dict(arrowstyle='->', color='red', alpha=0.7, lw=1),
                         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.9, edgecolor='red'),
                         fontsize=9, color='darkred', ha='center', fontweight='bold')
        root_annotations_D3.append(ann)
        offset_counter += 1
    
    for root in roots_D4:
        y_offset = -0.15 * 2 * y_range + 0.02 * 2 * y_range * (offset_counter % 3)
        ann = ax.annotate(f'D4: {root:.3f}', xy=(root, 0), 
                         xytext=(root, y_offset),
                         arrowprops=dict(arrowstyle='->', color='green', alpha=0.7, lw=1),
                         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.9, edgecolor='green'),
                         fontsize=9, color='darkgreen', ha='center', fontweight='bold')
        root_annotations_D4.append(ann)
        offset_counter += 1
    
    # Set annotation visibility based on current curve visibility
    for ann in root_annotations_D2:
        ann.set_visible(line_D2.get_visible())
    for ann in root_annotations_D3:
        ann.set_visible(line_D3.get_visible())
    for ann in root_annotations_D4:
        ann.set_visible(line_D4.get_visible())

def update_curve_visibility(label):
    """Update visibility of curves and corresponding roots"""
    # Toggle curve visibility based on checkbox state
    if label == 'Show D2':
        line_D2.set_visible(not line_D2.get_visible())
        scatter_D2.set_visible(line_D2.get_visible())
        # Update D2 root annotation visibility
        for ann in root_annotations_D2:
            ann.set_visible(line_D2.get_visible())
    elif label == 'Show D3':
        line_D3.set_visible(not line_D3.get_visible())
        scatter_D3.set_visible(line_D3.get_visible())
        # Update D3 root annotation visibility
        for ann in root_annotations_D3:
            ann.set_visible(line_D3.get_visible())
    elif label == 'Show D4':
        line_D4.set_visible(not line_D4.get_visible())
        scatter_D4.set_visible(line_D4.get_visible())
        # Update D4 root annotation visibility
        for ann in root_annotations_D4:
            ann.set_visible(line_D4.get_visible())
    
    # Redraw figure
    fig.canvas.draw_idle()

def update_plot():
    """Update the entire plot"""
    # Get current slider values
    sigma = sigma_slider.val
    s = s_slider.val
    x_min_val = xmin_slider.val
    x_max_val = xmax_slider.val
    
    # Ensure x_min < x_max
    if x_min_val >= x_max_val:
        x_min_val = x_max_val - 0.01
        xmin_slider.set_val(x_min_val)
    
    # Generate new x range
    m1_new = np.linspace(x_min_val, x_max_val, resolution)
    
    # Compute new values
    y_D2_new = np.array([D2(x, sigma, s) for x in m1_new])
    y_D3_new = np.array([D3(x, sigma, s) for x in m1_new])
    y_D4_new = np.array([D4(x, sigma, s) for x in m1_new])
    
    # Find new roots
    roots_D2_new = find_roots_in_range(D2, m1_new, y_D2_new)
    roots_D3_new = find_roots_in_range(D3, m1_new, y_D3_new)
    roots_D4_new = find_roots_in_range(D4, m1_new, y_D4_new)
    
    # Update curve data
    line_D2.set_xdata(m1_new)
    line_D2.set_ydata(y_D2_new)
    line_D3.set_xdata(m1_new)
    line_D3.set_ydata(y_D3_new)
    line_D4.set_xdata(m1_new)
    line_D4.set_ydata(y_D4_new)
    
    # Update root markers
    scatter_D2.set_offsets(np.column_stack([roots_D2_new, np.zeros(len(roots_D2_new))]))
    scatter_D3.set_offsets(np.column_stack([roots_D3_new, np.zeros(len(roots_D3_new))]))
    scatter_D4.set_offsets(np.column_stack([roots_D4_new, np.zeros(len(roots_D4_new))]))
    
    # Update root annotations
    update_root_annotations(roots_D2_new, roots_D3_new, roots_D4_new)
    
    # Remove old positive regions
    for patch in positive_patches:
        patch.remove()
    positive_patches.clear()
    
    # Find and plot new positive regions
    new_positive_regions = find_positive_intersection_regions(m1_new, y_D2_new, y_D3_new, y_D4_new)
    for region in new_positive_regions:
        patch = ax.axvspan(region[0], region[1], alpha=0.3, color='gold', zorder=0)
        positive_patches.append(patch)
    
    # Update legend for positive regions
    if positive_patches:
        positive_patches[0].set_label('All D>0')
    
    # Update axis limits
    ax.set_xlim(x_min_val, x_max_val)
    all_y_new = np.concatenate([y_D2_new, y_D3_new, y_D4_new])
    y_min_new, y_max_new = np.min(all_y_new), np.max(all_y_new)
    y_margin_new = 0.1 * (y_max_new - y_min_new) if y_max_new != y_min_new else 0.5
    ax.set_ylim(y_min_new - y_margin_new, y_max_new + y_margin_new)
    
    # Update title
    ax.set_title(f'D2-D4 Polynomials: σ={sigma:.2f}, s={s:.2f}, x∈[{x_min_val:.1f}, {x_max_val:.1f}]', 
                 fontsize=16, fontweight='bold')
    
    # Update info text
    update_info_text(sigma, s, x_min_val, x_max_val, 
                    roots_D2_new, roots_D3_new, roots_D4_new, new_positive_regions)
    
    # Update legend (keep lower left position)
    ax.legend(loc='lower left', fontsize=10, framealpha=0.9, bbox_to_anchor=(0.0, 0.0))
    
    # Redraw figure
    fig.canvas.draw_idle()

# Connect sliders to update function
sigma_slider.on_changed(lambda val: update_plot())
s_slider.on_changed(lambda val: update_plot())
xmin_slider.on_changed(lambda val: update_plot())
xmax_slider.on_changed(lambda val: update_plot())

# Connect checkboxes to visibility update function
check.on_clicked(update_curve_visibility)

# Initial update
update_info_text(sigma_initial, s_initial, x_min_initial, x_max_initial,
                roots_D2, roots_D3, roots_D4, positive_regions)

# Add shortcut instructions (lower left)
shortcut_text = ax.text(0.02, 0.02, 
                       'Controls:\nσ, s, x-range sliders\nCheckboxes: toggle curves',
                       transform=ax.transAxes, fontsize=8,
                       verticalalignment='bottom', horizontalalignment='left',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

plt.show()

# Print initial analysis
print("Initial analysis:")
print("=" * 60)
print(f"Parameters: σ={sigma_initial:.2f}, s={s_initial:.2f}")
print(f"x range: [{x_min_initial:.2f}, {x_max_initial:.2f}]")
print(f"Resolution: {resolution} points")
print(f"D2 roots ({len(roots_D2)}): {[f'{r:.4f}' for r in roots_D2]}")
print(f"D3 roots ({len(roots_D3)}): {[f'{r:.4f}' for r in roots_D3]}")
print(f"D4 roots ({len(roots_D4)}): {[f'{r:.4f}' for r in roots_D4]}")
print(f"All D>0 regions: {len(positive_regions)}")
if positive_regions:
    for i, region in enumerate(positive_regions):
        width = region[1] - region[0]
        print(f"  Region {i+1}: [{region[0]:.4f}, {region[1]:.4f}] (width={width:.4f})")
print("=" * 60)
print("Tips:")
print("1. All root positions are annotated with arrows, showing exact root locations")
print("2. Checkboxes at lower left control visibility of the three curves")
print("3. Sliders at lower right adjust parameters and x-axis range")
print("4. Info box at lower right shows current parameters, root positions, and positive regions")