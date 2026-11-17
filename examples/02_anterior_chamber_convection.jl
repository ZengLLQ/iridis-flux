"""
Example 02: 2D Sagittal Anterior Chamber Boussinesq Thermal Convection &
Intracameral Pigment Dispersion (Krukenberg Spindle) Simulation.

Demonstrates how the corneal-iris temperature differential drives an upward
thermal plume along the warm iris (37.0 C) and downward cooling current along
the corneal endothelium (34.0 C), causing inferior pigment sedimentation.
"""

using Printf
using IridisFlux

function main()
    println("====================================================================")
    println("EXAMPLE 02: ANTERIOR CHAMBER THERMAL CONVECTION & PIGMENT DISPERSION")
    println("====================================================================")

    healthy = healthy_emmetrope_phenotype()
    aacc = acute_angle_closure_phenotype()

    pde_open = solve_anterior_convection_2d(
        healthy.anatomy;
        nx = 33,
        ny = 49,
        t_cornea_C = 34.0,
        t_iris_C = 37.0,
        n_time_steps = 160,
        pigment_settling_mm_s = 0.025,
        iop_mmHg = 15.0,
    )

    pde_shallow = solve_anterior_convection_2d(
        aacc.anatomy;
        nx = 33,
        ny = 49,
        t_cornea_C = 34.0,
        t_iris_C = 37.0,
        n_time_steps = 160,
        pigment_settling_mm_s = 0.025,
        iop_mmHg = 46.0,
    )

    @printf("Open-Angle Chamber (3.10mm depth, 15 mmHg):\n")
    @printf("  Peak Convection Speed   : %.4f mm/s\n", pde_open.max_velocity_mm_s)
    @printf("  Corneal Nusselt Number  : %.3f\n", pde_open.mean_nusselt)
    @printf("  Inferior Pigment Ratio  : %.2f%%\n\n", 100.0 * pde_open.inferior_deposition_ratio)

    @printf("Shallow Angle-Closure Chamber (2.05mm depth, 46 mmHg):\n")
    @printf("  Peak Convection Speed   : %.4f mm/s\n", pde_shallow.max_velocity_mm_s)
    @printf("  Corneal Nusselt Number  : %.3f\n", pde_shallow.mean_nusselt)
    @printf("  Inferior Pigment Ratio  : %.2f%%\n\n", 100.0 * pde_shallow.inferior_deposition_ratio)

    println(render_anterior_chamber_map(pde_open))
end

main()
