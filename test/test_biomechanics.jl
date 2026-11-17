@testset "Lamina Cribrosa & Scleral Poroviscoelastic Biomechanics" begin
    anat = OcularAnatomy()
    bio = BiomechanicsParameters()

    @testset "Translaminar Pressure Gradient & Boundary Conditions" begin
        tlpg_norm = translaminar_pressure_gradient(15.0, 7.5, anat.lamina_thickness_mm)
        tlpg_high = translaminar_pressure_gradient(30.0, 7.5, anat.lamina_thickness_mm)
        tlpg_ntg = translaminar_pressure_gradient(15.0, 3.0, anat.lamina_thickness_mm)
        @test tlpg_high > tlpg_norm
        @test tlpg_ntg > tlpg_norm

        lam15 = evaluate_lamina_state(anat, bio, 15.0, 12.0)
        lam32 = evaluate_lamina_state(anat, bio, 32.0, 12.0)

        # Clamped at scleral canal rim (w(R) == 0), peak bowing at center (r = 0)
        @test isapprox(lam15.bowing_profile_um[end], 0.0; atol = 1e-10)
        @test lam15.peak_bowing_um == lam15.bowing_profile_um[1]
        @test lam32.peak_bowing_um > lam15.peak_bowing_um
        @test lam32.mean_axonal_strain > lam15.mean_axonal_strain
        @test lam32.axoplasmic_block_fraction > lam15.axoplasmic_block_fraction
        @test lam32.scleral_hoop_stress_kPa > lam15.scleral_hoop_stress_kPa
    end

    @testset "NTG Low CSF Pressure Elevation of Strain" begin
        ntg = normal_tension_glaucoma_phenotype()
        healthy = healthy_emmetrope_phenotype()
        # Same 16 mmHg IOP, but NTG has lower ICP and thinner/softer lamina
        lam_h = evaluate_lamina_state(healthy.anatomy, healthy.biomechanics, 16.0, 12.0)
        lam_n = evaluate_lamina_state(ntg.anatomy, ntg.biomechanics, 16.0, 12.0)
        @test lam_n.tlpg_mmHg_mm > lam_h.tlpg_mmHg_mm
        @test lam_n.peak_bowing_um > lam_h.peak_bowing_um
        @test lam_n.mean_axonal_strain > lam_h.mean_axonal_strain
    end
end
