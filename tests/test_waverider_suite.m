function results = test_waverider_suite()
%TEST_WAVERIDER_SUITE Automated unit and manifold watertightness test suite.
%   Verifies Taylor-Maccoll conical flow accuracy, oblique shock relations,
%   all four waverider generation methodologies ('conical', 'osculating_cone',
%   'caret', 'cst'), both upper surface modes ('freestream', 'canopy'), and
%   STL/OBJ mesh export integrity.

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootDir);

    fprintf('================================================================\n');
    fprintf('  RUNNING HYPERSONIC WAVERIDER AUTOMATED TEST SUITE\n');
    fprintf('================================================================\n');

    % 1. Test Oblique Shock & Taylor-Maccoll Conical Flow Solver
    shock = waverider.flow.obliqueShock(6.0, 15.0, 1.4);
    assert(abs(shock.deltaDeg - 7.18494) < 1e-4, 'Oblique shock deflection mismatch.');
    assert(abs(shock.p2_p1 - 2.64680) < 1e-4, 'Oblique shock pressure ratio mismatch.');

    tm = waverider.flow.taylorMaccoll(6.0, 15.0, 1.4);
    assert(abs(tm.thetaConeDeg - 10.73711) < 2e-3, 'Taylor-Maccoll cone angle mismatch.');
    assert(abs(tm.utheta(1)) < 1e-9, 'Normal velocity at cone surface must be zero.');
    fprintf('[PASS] Flowfield solvers: M=6, beta=15 deg -> delta_wedge=%.4f deg, theta_cone=%.4f deg\n', ...
        shock.deltaDeg, tm.thetaConeDeg);

    % Verify sub-Mach shock angle guard
    caughtSubMach = false;
    try
        waverider.WaveriderBody('Mach', 3.0, 'ShockAngleDeg', 15.0); % asin(1/3) = 19.47 deg
    catch ME
        caughtSubMach = contains(ME.identifier, 'ShockBelowMachAngle');
    end
    assert(caughtSubMach, 'Expected ShockBelowMachAngle error when beta <= asin(1/M).');

    % 2. Test all 4 Waverider Presets + Both Upper Surface Types
    presets = {'conical', 'osculating_cone', 'caret', 'cst'};
    upperModes = {'freestream', 'canopy'};

    results = struct();
    idx = 0;

    for p = 1:numel(presets)
        for u = 1:numel(upperModes)
            idx = idx + 1;
            cfg = waverider.WaveriderConfig.fromPreset(presets{p});
            cfg.UpperSurfaceType = upperModes{u};
            if strcmp(upperModes{u}, 'canopy') && cfg.CanopyHeightRatio == 0
                cfg.CanopyHeightRatio = 0.05;
            end

            body = waverider.WaveriderBody(cfg);
            rep = body.Report;

            assert(rep.IsWatertight, ...
                sprintf('Mesh failed watertight verification for %s (%s).', presets{p}, upperModes{u}));
            assert(rep.NumBoundaryEdges == 0, 'Mesh has open free boundary edges.');
            assert(rep.NumNonManifoldEdges == 0, 'Mesh has non-manifold edges.');
            assert(rep.NumConflictingHalfEdges == 0, 'Mesh has inconsistent face normal orientations.');
            assert(rep.EulerCharacteristic == 2, 'Closed genus-0 mesh must have Euler characteristic = 2.');
            assert(rep.EnclosedVolume > 0.1, 'Enclosed volume must be strictly positive.');
            assert(body.AeroMetrics.CL > 0 && body.AeroMetrics.LD_inviscid > 1.5, ...
                'Aerodynamic L/D out of expected physical bounds.');

            fprintf('[PASS] %-16s | Upper: %-10s | Faces: %5d | Vol: %6.3f m^3 | eta_V: %.4f | L/D: %.2f | Euler: %d\n', ...
                presets{p}, upperModes{u}, rep.NumFaces, rep.EnclosedVolume, ...
                rep.VolumetricEfficiency, body.AeroMetrics.LD_inviscid, rep.EulerCharacteristic);

            results(idx).Preset = presets{p};
            results(idx).UpperMode = upperModes{u};
            results(idx).Report = rep;
            results(idx).Aero = body.AeroMetrics;
        end
    end

    % 3. Test Conformal Aft Frustum & Conical Rocket Booster Mounts
    blendLaws = {'linear', 'cubic', 'quintic', 'conic'};
    for b = 1:numel(blendLaws)
        idx = idx + 1;
        cfgF = waverider.WaveriderConfig.fromPreset('booster_mounted');
        cfgF.FrustumBlendLaw = blendLaws{b};
        cfgF.EnableRocketStage = (b > 1); % Test frustum-only (b==1) and frustum+rocket (b>1)

        bodyF = waverider.WaveriderBody(cfgF);
        repF = bodyF.Report;

        assert(repF.IsWatertight, sprintf('Frustum mount (%s) failed watertight check.', blendLaws{b}));
        assert(repF.NumBoundaryEdges == 0, 'Frustum mesh has open boundary edges.');
        assert(repF.NumNonManifoldEdges == 0, 'Frustum mesh has non-manifold edges.');
        assert(repF.NumConflictingHalfEdges == 0, 'Frustum mesh has conflicting half-edges.');
        assert(repF.EulerCharacteristic == 2, 'Frustum mesh must have Euler characteristic = 2.');

        fprintf('[PASS] booster_mounted  | Blend: %-8s | Faces: %5d | Vol: %6.3f m^3 | Rocket: %d | Euler: %d\n', ...
            blendLaws{b}, repF.NumFaces, repF.EnclosedVolume, double(cfgF.EnableRocketStage), repF.EulerCharacteristic);

        results(idx).Preset = 'booster_mounted';
        results(idx).UpperMode = blendLaws{b};
        results(idx).Report = repF;
        results(idx).Aero = bodyF.AeroMetrics;
    end

    % 4. Test STL and OBJ Export and Read-Back Verification
    outDir = fullfile(rootDir, 'output');
    bodyOC = waverider.WaveriderBody('osculating_cone');
    stlFile = bodyOC.exportSTL(fullfile(outDir, 'test_osculating_cone.stl'));
    objFile = bodyOC.exportOBJ(fullfile(outDir, 'test_osculating_cone.obj'));

    assert(isfile(stlFile), 'STL file was not created.');
    assert(isfile(objFile), 'OBJ file was not created.');

    trRead = stlread(stlFile);
    fbRead = freeBoundary(trRead);
    assert(isempty(fbRead), 'Re-imported STL triangulation is not watertight!');
    objText = fileread(objFile);
    nObjV = numel(regexp(objText, '(?m)^v\s+'));
    nObjF = numel(regexp(objText, '(?m)^f\s+'));
    assert(nObjV == bodyOC.Report.NumVertices && nObjF == bodyOC.Report.NumFaces, ...
        'OBJ vertex/face count does not match mesh report.');
    pol0 = bodyOC.computeAeroPolar(0.0);
    assert(abs(pol0.CL - bodyOC.AeroMetrics.CL) < 1e-13, 'AeroPolar(0) must match AeroMetrics.CL.');
    fprintf('[PASS] STL & OBJ export verified (STL freeBoundary = %d, OBJ V=%d, F=%d)\n', ...
        size(fbRead, 1), nObjV, nObjF);

    % 5. Test Range-Angle-Intensity (RTI) & Waveform Bandwidth Resolution
    rti200  = waverider.rcs.computeRangeProfile(bodyOC, 10.0, 200.0,  [0, 45, 90], -22.0);
    rti1000 = waverider.rcs.computeRangeProfile(bodyOC, 10.0, 1000.0, [0, 45, 90], -22.0);
    assert(abs(rti200.deltaR_m - 5.0 * rti1000.deltaR_m) < 1e-6, 'Range resolution must scale inversely with bandwidth B.');
    assert(all(isfinite(rti1000.RTI_total_dBsm(:))), 'RTI matrix contains non-finite values.');
    fprintf('[PASS] Range-Angle-Intensity verified (B=1000 MHz -> DeltaR = %.2f cm, B=200 MHz -> DeltaR = %.2f cm)\n', ...
        rti1000.deltaR_m * 100, rti200.deltaR_m * 100);

    % 6. Test Conventional Reentry Vehicle (RV) Presets & Polarimetric RCS
    rvPresets = {'rv_monoconic', 'rv_biconic', 'rv_triconic_flare', 'rv_sharp_conic'};
    for p = 1:numel(rvPresets)
        idx = idx + 1;
        cfgRV = waverider.WaveriderConfig.fromPreset(rvPresets{p});
        bodyRV = waverider.WaveriderBody(cfgRV);
        repRV = bodyRV.Report;

        assert(repRV.IsWatertight, sprintf('RV preset %s failed watertight check.', rvPresets{p}));
        assert(repRV.NumBoundaryEdges == 0, 'RV mesh has open boundary edges.');
        assert(repRV.NumNonManifoldEdges == 0, 'RV mesh has non-manifold edges.');
        assert(repRV.NumConflictingHalfEdges == 0, 'RV mesh has conflicting half-edges.');
        assert(repRV.EulerCharacteristic == 2, 'RV mesh must have Euler characteristic = 2.');
        assert(abs(bodyRV.AeroMetrics.CL) < 1e-13 && abs(bodyRV.AeroMetrics.CY) < 1e-13, ...
            'Axisymmetric RV at alpha=0 must have zero lift and zero side force.');

        rcsRV = bodyRV.computeRCS(10.0, [0; 90; 170], 0.0, 'Polarization', 'OP');
        assert(all(isfinite(rcsRV.SigmaTotal_dBsm(:))), 'RV polarimetric RCS returned non-finite values.');

        fprintf('[PASS] %-18s | Nose: %-9s | Aft: %-7s | Faces: %5d | Vol: %6.3f m^3 | CDw: %.4f | Euler: %d\n', ...
            rvPresets{p}, cfgRV.RVNoseType, cfgRV.RVAftType, repRV.NumFaces, ...
            repRV.EnclosedVolume, bodyRV.AeroMetrics.CD_wave, repRV.EulerCharacteristic);

        results(idx).Preset = rvPresets{p};
        results(idx).UpperMode = cfgRV.RVAftType;
        results(idx).Report = repRV;
        results(idx).Aero = bodyRV.AeroMetrics;
    end

    fprintf('================================================================\n');
    fprintf('  ALL %d WAVERIDER & RV TEST CONFIGURATIONS PASSED 100%%!\n', idx);
    fprintf('================================================================\n');
end
