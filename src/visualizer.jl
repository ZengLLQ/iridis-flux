"""
Pure 7-bit ASCII terminal visualizers for IridisFlux.

Strictly uses 7-bit ASCII characters (codes 0-127: `+`, `-`, `|`, `*`, `#`,
`.`, `:`, `^`, `v`, `o`, `@`, `=`) so Windows cp1252 and headless CI consoles
never raise encoding errors.

Provides:
- `render_iop_opp_chart`: Dual time-series trajectory of IOP and OPP.
- `render_anterior_chamber_map`: 2D sagittal thermal convection & pigment density map.
- `render_lamina_profile`: Radial lamina cribrosa posterior bowing profile.
- `render_visual_field_24_2`: Humphrey 24-2 54-point sensitivity grid (dB & grayscale).
- `render_clinical_card`: Structured clinical diagnostic & progression summary box.
"""

using Printf

function _ascii_border(width::Int)
    return "+" * repeat("-", width - 2) * "+"
end

function _ascii_row(text::String, width::Int)
    inner = width - 4
    clipped = length(text) > inner ? text[1:inner] : text * repeat(" ", inner - length(text))
    return "| " * clipped * " |"
end

"""
    render_iop_opp_chart(res; width=68, height=14)

Render a dual-trace 7-bit ASCII time-series chart showing Intraocular Pressure
(`*` = IOP mmHg) and Ocular Perfusion Pressure (`o` = OPP mmHg) across the
simulation horizon.
"""
function render_iop_opp_chart(res::DiurnalSimulationResult; width::Int = 68, height::Int = 14)
    lines = String[]
    push!(lines, _ascii_border(width))
    title = @sprintf("DIURNAL IOP (*) & OPP (o) TRAJECTORY | %s | %s", res.phenotype_name, res.regimen_name)
    push!(lines, _ascii_row(title, width))
    push!(lines, _ascii_border(width))

    plot_w = width - 14
    all_vals = vcat(res.iop_mmHg, res.opp_mmHg)
    y_min = floor(minimum(all_vals) / 5.0) * 5.0
    y_max = ceil(maximum(all_vals) / 5.0) * 5.0
    if y_max <= y_min + 5.0
        y_max = y_min + 10.0
    end

    canvas = fill(' ', height, plot_w)
    # Reference 21 mmHg ocular hypertension line
    if y_min <= 21.0 <= y_max
        row_21 = clamp(height - round(Int, (21.0 - y_min) / (y_max - y_min) * (height - 1)), 1, height)
        for c in 1:plot_w
            canvas[row_21, c] = '.'
        end
    end

    n_pts = length(res.time_h)
    for col in 1:plot_w
        idx = clamp(1 + round(Int, (col - 1) / max(plot_w - 1, 1) * (n_pts - 1)), 1, n_pts)
        r_opp = clamp(height - round(Int, (res.opp_mmHg[idx] - y_min) / (y_max - y_min) * (height - 1)), 1, height)
        r_iop = clamp(height - round(Int, (res.iop_mmHg[idx] - y_min) / (y_max - y_min) * (height - 1)), 1, height)
        canvas[r_opp, col] = 'o'
        canvas[r_iop, col] = (r_iop == r_opp) ? '@' : '*'
    end

    for r in 1:height
        y_val = y_max - (r - 1) / max(height - 1, 1) * (y_max - y_min)
        label = @sprintf("%5.1f |", y_val)
        row_str = label * String(canvas[r, :])
        push!(lines, _ascii_row(row_str, width))
    end
    axis = "      +" * repeat("-", plot_w)
    push!(lines, _ascii_row(axis, width))
    t0 = res.time_h[1]
    tf = res.time_h[end]
    tmid = 0.5 * (t0 + tf)
    t_lbl = @sprintf("      %-8.1fh%s%-8.1fh%s%6.1fh", t0, repeat(" ", div(plot_w, 2) - 10), tmid, repeat(" ", plot_w - div(plot_w, 2) - 14), tf)
    push!(lines, _ascii_row(t_lbl, width))
    legend = "Legend: * IOP (mmHg)  o OPP (mmHg)  . 21 mmHg threshold"
    push!(lines, _ascii_row(legend, width))
    push!(lines, _ascii_border(width))
    return join(lines, "\n")
end

"""
    render_anterior_chamber_map(pde; width=68)

Render a 2D sagittal heatmap of anterior chamber aqueous thermal convection
velocity & intracameral pigment deposition (`x` = Cornea -> Iris, `y` = Inferior -> Superior).
"""
function render_anterior_chamber_map(pde::AnteriorChamberPDEState; width::Int = 68)
    lines = String[]
    push!(lines, _ascii_border(width))
    push!(lines, _ascii_row("ANTERIOR CHAMBER BOUSSINESQ CONVECTION & PIGMENT DISPERSION", width))
    sub = @sprintf("Max Vel: %.3f mm/s | Nu: %.2f | Inferior Deposition: %.1f%%", pde.max_velocity_mm_s, pde.mean_nusselt, 100.0 * pde.inferior_deposition_ratio)
    push!(lines, _ascii_row(sub, width))
    push!(lines, _ascii_border(width))

    palette = " .:-=+*#%@"
    n_pal = length(palette)
    rows_out = 16
    cols_out = 22

    c_max = max(maximum(pde.tracer_conc), 1e-12)
    t_min = minimum(pde.temperature_C)
    t_max = max(maximum(pde.temperature_C), t_min + 1e-6)

    for r in 1:rows_out
        # Top of plot = Superior (y = H), Bottom = Inferior (y = 0)
        j = clamp(round(Int, (1.0 - (r - 1) / (rows_out - 1)) * (pde.ny - 1)) + 1, 1, pde.ny)
        buf_c = IOBuffer()
        buf_t = IOBuffer()
        for c in 1:cols_out
            i = clamp(round(Int, (c - 1) / (cols_out - 1) * (pde.nx - 1)) + 1, 1, pde.nx)
            norm_c = clamp(pde.tracer_conc[i, j] / c_max, 0.0, 1.0)
            idx_c = clamp(round(Int, norm_c * (n_pal - 1)) + 1, 1, n_pal)
            print(buf_c, palette[idx_c])

            norm_t = clamp((pde.temperature_C[i, j] - t_min) / (t_max - t_min), 0.0, 1.0)
            idx_t = clamp(round(Int, norm_t * (n_pal - 1)) + 1, 1, n_pal)
            print(buf_t, palette[idx_t])
        end
        c_str = String(take!(buf_c))
        t_str = String(take!(buf_t))
        line_s = @sprintf("y=%4.1f |%s|   |%s|", pde.y_mm[j], t_str, c_str)
        push!(lines, _ascii_row(line_s, width))
    end
    push!(lines, _ascii_row("       Cornea -> Iris (Temp)     Cornea -> Iris (Pigment)", width))
    push!(lines, _ascii_border(width))
    return join(lines, "\n")
end

"""
    render_lamina_profile(lam; width=68)

Render radial lamina cribrosa posterior bowing profile `w(r)` (um) from
optic disc center (`r = 0`) to scleral canal rim (`r = R_LC`).
"""
function render_lamina_profile(lam::LaminaState; width::Int = 68)
    lines = String[]
    push!(lines, _ascii_border(width))
    hdr = @sprintf("LAMINA CRIBROSA BIOMECHANICS | TLPG=%.1f mmHg/mm | Peak=%.1f um", lam.tlpg_mmHg_mm, lam.peak_bowing_um)
    push!(lines, _ascii_row(hdr, width))
    push!(lines, _ascii_border(width))

    bar_max = width - 28
    n_pts = length(lam.radial_coords_mm)
    w_peak = max(maximum(abs.(lam.bowing_profile_um)), 1e-3)
    step = max(1, div(n_pts, 9))
    for i in 1:step:n_pts
        w_val = lam.bowing_profile_um[i]
        len = clamp(round(Int, abs(w_val) / w_peak * bar_max), 0, bar_max)
        bar = repeat("#", len) * repeat(".", bar_max - len)
        row = @sprintf("r=%.2fmm [%s] %5.1fum", lam.radial_coords_mm[i], bar, w_val)
        push!(lines, _ascii_row(row, width))
    end
    ftr = @sprintf("Axonal Strain: %.4f | BDNF Block: %.1f%% | Hoop: %.1f kPa", lam.mean_axonal_strain, 100.0 * lam.axoplasmic_block_fraction, lam.scleral_hoop_stress_kPa)
    push!(lines, _ascii_row(ftr, width))
    push!(lines, _ascii_border(width))
    return join(lines, "\n")
end

"""
    render_visual_field_24_2(vf; width=68)

Render the 54-point Humphrey 24-2 visual field sensitivity grid (numeric dB
on left, defect density symbols on right) with global perimetric indices.
"""
function render_visual_field_24_2(vf::VisualFieldMap; width::Int = 68)
    lines = String[]
    push!(lines, _ascii_border(width))
    hdr = @sprintf("HUMPHREY 24-2 VISUAL FIELD | MD=%+.2f dB | PSD=%.2f | VFI=%.1f%%", vf.mean_deviation_dB, vf.pattern_std_dB, vf.visual_field_index_pct)
    push!(lines, _ascii_row(hdr, width))
    push!(lines, _ascii_border(width))

    for r in 1:8
        buf_num = IOBuffer()
        buf_sym = IOBuffer()
        for c in 1:9
            if vf.valid_mask[r, c]
                db = vf.grid_dB[r, c]
                @printf(buf_num, "%3d", round(Int, db))
                sym = if db >= 28.0
                    " ."
                elseif db >= 22.0
                    " :"
                elseif db >= 15.0
                    " +"
                elseif db >= 8.0
                    " *"
                else
                    " #"
                end
                print(buf_sym, sym)
            else
                print(buf_num, "   ")
                print(buf_sym, (r in (4, 5) && c == 8) ? " B" : "  ")
            end
        end
        row = String(take!(buf_num)) * "  |  " * String(take!(buf_sym))
        push!(lines, _ascii_row(row, width))
    end
    push!(lines, _ascii_row("Defect map: . >=28dB  : 22-27  + 15-21  * 8-14  # <8dB  B=BlindSpot", width))
    push!(lines, _ascii_border(width))
    return join(lines, "\n")
end

"""
    render_clinical_card(patient, diurnal_res, prog_res, risk; width=68)

Render a compact multi-metric 7-bit ASCII diagnostic summary card combining
hemodynamic, biomechanical, pharmacological, and 5-year visual field outcomes.
"""
function render_clinical_card(
    patient::PatientPhenotype,
    diurnal_res::DiurnalSimulationResult,
    prog_res::LongTermProgressionResult,
    risk;
    width::Int = 68,
)
    lines = String[]
    push!(lines, _ascii_border(width))
    push!(lines, _ascii_row("IRIDISFLUX CLINICAL DIAGNOSTIC & PROGNOSIS REPORT", width))
    push!(lines, _ascii_border(width))
    push!(lines, _ascii_row("Patient Phenotype : " * patient.clinical_label, width))
    push!(lines, _ascii_row("Therapy Regimen   : " * diurnal_res.regimen_name, width))
    push!(lines, _ascii_row("Risk Stage        : " * risk.stage, width))
    push!(lines, _ascii_row("Primary Driver    : " * risk.dominant_driver, width))
    push!(lines, _ascii_border(width))
    push!(lines, _ascii_row(@sprintf("IOP (Mean / Peak / Swing): %.2f / %.2f / %.2f mmHg", diurnal_res.mean_iop, diurnal_res.peak_iop, diurnal_res.iop_fluctuation), width))
    push!(lines, _ascii_row(@sprintf("Minimum Nocturnal OPP    : %.2f mmHg", diurnal_res.min_opp), width))
    push!(lines, _ascii_row(@sprintf("Peak Lamina Bowing       : %.2f um (Strain: %.4f)", maximum(diurnal_res.lamina_bowing_um), diurnal_res.peak_strain), width))
    push!(lines, _ascii_row(@sprintf("Mean ONH Hypoxic Stress  : %.3f (Min PO2: %.1f mmHg)", diurnal_res.mean_hypoxic_stress, minimum(diurnal_res.min_po2_mmHg)), width))
    push!(lines, _ascii_row(@sprintf("5-Year RNFL Thickness    : %.1f um (Start: %.1f um)", prog_res.mean_rnfl_um[end], prog_res.mean_rnfl_um[1]), width))
    push!(lines, _ascii_row(@sprintf("5-Year 24-2 MD / VFI     : %+.2f dB (%.2f dB/yr) | VFI %.1f%%", risk.md_5yr_dB, risk.md_rate_dB_yr, risk.vfi_5yr_pct), width))
    push!(lines, _ascii_border(width))
    return join(lines, "\n")
end
