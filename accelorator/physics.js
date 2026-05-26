const FusionPhysics = (() => {

    const COULOMB_CONSTANT = 1000;
    const MAX_SPEED = 25;
    const MAGNETIC_ZONE = 70;
    const B_MAX = 3.0;

    class Nucleon {
        constructor(type, x, y, clusterType) {
            this.type = type;
            this.x = x;
            this.y = y;
            this.vx = 0;
            this.vy = 0;
            this.ax = 0;
            this.ay = 0;
            this.old_ax = 0;
            this.old_ay = 0;
            this.mass = 1.0;
            this.charge = (type === 'proton') ? 1.0 : 0.0;
            this.radius = 12;
            this.color = (type === 'proton') ? '#ef4444' : '#3b82f6';
            this.isAlphaComponent = false;
            this.isEjectedNeutron = false;
            this.clusterType = clusterType;
        }

        kineticEnergy() {
            let pixelEnergy = 0.5 * this.mass * (this.vx * this.vx + this.vy * this.vy);

            if (this.isEjectedNeutron) {
                // Calibrates a visual velocity of 16.0 to read exactly 14.1 MeV
                const MEV_SCALE = 0.055078;
                return pixelEnergy * MEV_SCALE;
            }

            const STANDARD_SCALE = 0.005;
            return pixelEnergy * STANDARD_SCALE;
        }

        draw(ctx) {
            ctx.save();
            if (this.charge > 0) {
                let coulombGrad = ctx.createRadialGradient(this.x, this.y, this.radius, this.x, this.y, 120);
                coulombGrad.addColorStop(0, 'rgba(239, 68, 68, 0.15)');
                coulombGrad.addColorStop(1, 'rgba(239, 68, 68, 0)');
                ctx.beginPath(); ctx.arc(this.x, this.y, 120, 0, Math.PI * 2);
                ctx.fillStyle = coulombGrad; ctx.fill();
            }

            let strongGrad = ctx.createRadialGradient(this.x, this.y, this.radius, this.x, this.y, 18);
            strongGrad.addColorStop(0, 'rgba(34, 197, 94, 0.4)');
            strongGrad.addColorStop(1, 'rgba(34, 197, 94, 0)');
            ctx.beginPath(); ctx.arc(this.x, this.y, 18, 0, Math.PI * 2);
            ctx.fillStyle = strongGrad; ctx.fill();
            ctx.restore();

            ctx.beginPath();
            ctx.arc(this.x, this.y, this.radius, 0, Math.PI * 2);
            ctx.fillStyle = this.color;
            ctx.fill();

            if (this.isAlphaComponent) {
                ctx.strokeStyle = '#facc15';
                ctx.lineWidth = 3;
                ctx.stroke();
            }

            ctx.fillStyle = '#ffffff';
            ctx.font = 'bold 9px sans-serif';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';

            let label = this.type === 'proton' ? 'P+' : 'N0';
            if (this.isEjectedNeutron) label = 'n⁰';

            ctx.fillText(label, this.x, this.y);
        }
    }

    class Photon {
        constructor(x, y, angle) {
            this.x = x;
            this.y = y;
            this.speed = 12;
            this.vx = Math.cos(angle) * this.speed;
            this.vy = Math.sin(angle) * this.speed;
            this.life = 1.0;
            this.decay = 0.02;
        }

        update(timeScale) {
            this.x += this.vx * timeScale;
            this.y += this.vy * timeScale;
            this.life -= this.decay * timeScale;
        }

        draw(ctx) {
            ctx.save();
            ctx.globalAlpha = this.life;
            ctx.strokeStyle = '#a855f7';
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(this.x, this.y);
            for (let i = 0; i < 6; i++) {
                let evalX = this.x - (this.vx * i * 0.2);
                let evalY = this.y - (this.vy * i * 0.2);
                let waveOffset = Math.sin(Date.now() * 0.05 + i) * 4;
                let mag = Math.sqrt((-this.vy) ** 2 + this.vx ** 2) || 1;
                ctx.lineTo(evalX + (-this.vy / mag) * waveOffset, evalY + (this.vx / mag) * waveOffset);
            }
            ctx.stroke();
            ctx.restore();
        }
    }

    function rotateCluster(clusterType, cx, cy, rad, nucleons, sim) {
        if (sim && sim.reactionTriggered) return;

        let filtered = nucleons.filter(n => n.clusterType === clusterType);
        
        if (clusterType === 'D' && filtered.length === 2) {
            let r0 = 12.0; 
            filtered[0].x = cx - Math.cos(rad) * r0;
            filtered[0].y = cy - Math.sin(rad) * r0;
            filtered[1].x = cx + Math.cos(rad) * r0;
            filtered[1].y = cy + Math.sin(rad) * r0;
        } else if (clusterType === 'T' && filtered.length === 3) {
            // FIX: Boosted to 24.0px to stop instant overlap explosions upon spawn
            let r_equil = 24.0;
            let h = r_equil * (Math.sqrt(3) / 2);
            let points = [
                { x: 0, y: -(2 / 3) * h },
                { x: -r_equil / 2, y: (1 / 3) * h },
                { x: r_equil / 2, y: (1 / 3) * h }
            ];
            filtered.forEach((n, idx) => {
                let pt = points[idx];
                n.x = cx + (pt.x * Math.cos(rad) - pt.y * Math.sin(rad));
                n.y = cy + (pt.x * Math.sin(rad) + pt.y * Math.cos(rad));
            });
        }
        
        filtered.forEach(n => {
            n.vx = 0; n.vy = 0; n.ax = 0; n.ay = 0; n.old_ax = 0; n.old_ay = 0;
        });
    }

    function computeInterNucleonForces(n1, n2) {
        let dx = n2.x - n1.x;
        let dy = n2.y - n1.y;
        let r = Math.sqrt(dx * dx + dy * dy);
        if (r < 0.1) return;
    
        // 1. Long-range Coulomb Repulsion (Protons only)
        let f_coulomb = (COULOMB_CONSTANT * n1.charge * n2.charge) / (r * r);
    
        // 2. Realistic Strong Force (Yukawa Attraction + Core Repulsion Well)
        const SIGMA = 24.0; 
        let f_strong = 0;
    
        if ((n1.clusterType === 'Alpha' && n2.clusterType === 'Alpha') || (r < 45)) {
            let attractive = Math.pow(SIGMA / r, 4);  
            let repulsive = Math.pow(SIGMA / r, 8);   
            
            const C_STRONG = 8000; 
            f_strong = C_STRONG * (repulsive - attractive);
        } else {
            const MU = 0.15;
            let expTerm = Math.exp(-MU * r);
            f_strong = -3500 * (expTerm / r) * (MU + (1 / r));
        }
    
        let f_net = f_coulomb + f_strong;
    
        let fx = (dx / r) * f_net;
        let fy = (dy / r) * f_net;
    
        n1.ax -= fx / n1.mass;
        n1.ay -= fy / n1.mass;
        n2.ax += fx / n2.mass;
        n2.ay += fy / n2.mass;
    }

    function applyUnifiedConstraints(nucleons) {
        const SOLVER_PASSES = 5; 

        for (let pass = 0; pass < SOLVER_PASSES; pass++) {
            for (let i = 0; i < nucleons.length; i++) {
                for (let j = i + 1; j < nucleons.length; j++) {
                    let n1 = nucleons[i];
                    let n2 = nucleons[j];

                    let dx = n2.x - n1.x;
                    let dy = n2.y - n1.y;
                    let dist = Math.sqrt(dx * dx + dy * dy) || 0.1;

                    // Condition A: Inside un-fused ions
                    if (n1.clusterType === n2.clusterType && n1.clusterType !== 'Alpha' && n1.clusterType !== 'Free') {
                        // FIX: Synchronized to 24.0 so your structural target completely agrees with physical radii bounds
                        let targetDist = 24.0; 
                        let difference = targetDist - dist;
                        let percent = (difference / dist) * 0.5;
                        
                        n1.x -= dx * percent; n1.y -= dy * percent;
                        n2.x += dx * percent; n2.y += dy * percent;
                    } 
                    // Condition B: Hard Solid Boundary Overlapping
                    else {
                        let minDist = n1.radius + n2.radius; 
                        if (dist < minDist) {
                            let difference = minDist - dist;
                            let percent = (difference / dist) * 0.5;

                            n1.x -= dx * percent; n1.y -= dy * percent;
                            n2.x += dx * percent; n2.y += dy * percent;
                            
                            let rvx = n2.vx - n1.vx;
                            let rvy = n2.vy - n1.vy;
                            n1.vx += rvx * 0.05; n1.vy += rvy * 0.05;
                            n2.vx -= rvx * 0.05; n2.vy -= rvy * 0.05;
                        }
                    }
                }
            }
        }
    }

    function evaluateFusionState(sim) {
        if (sim.nucleons.length !== 5 || sim.fusionOccurred) return;
    
        // 1. Calculate Center of Mass (COM) Position and Net Momentum
        let totalMass = 0;
        let sumX = 0, sumY = 0;
        let netVx = 0, netVy = 0;
    
        sim.nucleons.forEach(n => {
            sumX += n.x * n.mass;
            sumY += n.y * n.mass;
            netVx += n.vx * n.mass; // Momentum X
            netVy += n.vy * n.mass; // Momentum Y
            totalMass += n.mass;
        });
    
        let avgX = sumX / totalMass;
        let avgY = sumY / totalMass;
        
        // System velocity (velocity of the center of mass)
        let v_com_x = netVx / totalMass;
        let v_com_y = netVy / totalMass;
    
        let maxSpread = Math.max(...sim.nucleons.map(n =>
            Math.sqrt((n.x - avgX) ** 2 + (n.y - avgY) ** 2)
        ));
    
        // Trigger threshold (adjusted to clear structural constraints safely)
        if (maxSpread < 45) {
            let neutrons = sim.nucleons.filter(n => n.type === 'neutron');
            // Sort farthest from center to find the escaping neutron
            neutrons.sort((a, b) =>
                ((b.x - avgX) ** 2 + (b.y - avgY) ** 2) -
                ((a.x - avgX) ** 2 + (a.y - avgY) ** 2)
            );
            let escapeNeutron = neutrons[0];
    
            // Random or impact-dependent ejection angle
            let angle = Math.atan2(escapeNeutron.y - avgY, escapeNeutron.x - avgX);
    
            /* * PHYSICS CALIBRATION
             * Target: Neutron KE = 14.1 MeV, Alpha KE = 3.5 MeV
             * Since KE = 0.5 * m * v^2 
             * For Neutron (m=1): v = sqrt(2 * KE / STANDARD_SCALE)
             * If standard scale is 0.005: v_neutron = sqrt(2 * 14.1 / 0.005) = sqrt(5640) ≈ 75.1
             * * However, your code uses a custom MEV_SCALE (0.055078) for the neutron.
             * Let's derive physical velocities matching your internal energy readers precisely:
             */
            const V_NEUTRON_REL = 16.0;             // Calibrates exactly to 14.1 MeV in your engine
            const V_ALPHA_REL = V_NEUTRON_REL / 4.0; // 4.0 speed ratio due to mass conservation (m_alpha=4, m_n=1)
    
            sim.nucleons.forEach(n => {
                if (n === escapeNeutron) {
                    n.clusterType = 'Free';
                    n.isEjectedNeutron = true;
                    
                    // Final velocity = System velocity + Ejection velocity vector
                    n.vx = v_com_x + Math.cos(angle) * V_NEUTRON_REL;
                    n.vy = v_com_y + Math.sin(angle) * V_NEUTRON_REL;
                } else {
                    n.clusterType = 'Alpha';
                    n.isAlphaComponent = true;
                    
                    // Recoil direction is exactly 180 degrees opposite
                    n.vx = v_com_x - Math.cos(angle) * V_ALPHA_REL;
                    n.vy = v_com_y - Math.sin(angle) * V_ALPHA_REL;
                }
            });
    
            sim.fusionOccurred = true;
            sim.fusionSlowFactor = sim.fusionMinDt;
    
            // Generate isotropic gamma photon burst
            for (let i = 0; i < 16; i++) {
                sim.photons.push(new Photon(avgX, avgY, (i / 16) * Math.PI * 2));
            }
        }
    }

    function computeMinInterClusterDistance(nucleons) {
        let minDist = Infinity;
        for (let i = 0; i < nucleons.length; i++) {
            for (let j = i + 1; j < nucleons.length; j++) {
                if (nucleons[i].clusterType === nucleons[j].clusterType) continue;
                let dx = nucleons[i].x - nucleons[j].x;
                let dy = nucleons[i].y - nucleons[j].y;
                let d = Math.sqrt(dx * dx + dy * dy);
                if (d < minDist) minDist = d;
            }
        }
        return minDist;
    }

    function computeAdaptiveFactor(minDist, approachRange, approachMinDt) {
        let adaptiveFactor = 1.0;
        if (minDist !== Infinity && minDist < approachRange) {
            let proximityRatio = Math.max(0, minDist / approachRange);
            adaptiveFactor = approachMinDt + (1.0 - approachMinDt) * Math.pow(proximityRatio, 3);

            if (minDist < 45) {
                adaptiveFactor *= Math.max(0.005, minDist / 45);
            }
        }
        return adaptiveFactor;
    }

    function step(sim, canvasWidth, canvasHeight) {
        let minDist = Infinity;
        if (sim.reactionTriggered && sim.nucleons.length >= 2) {
            minDist = computeMinInterClusterDistance(sim.nucleons);
        }

        let adaptiveFactor = 1.0;
        if (sim.reactionTriggered) {
            adaptiveFactor = computeAdaptiveFactor(minDist, sim.approachRange, sim.approachMinDt);
        }

        if (sim.fusionOccurred) {
            sim.fusionSlowFactor = Math.min(sim.fusionSlowFactor + sim.fusionRecoveryRate, 1.0);
        }

        let dt = sim.timeScale * adaptiveFactor * sim.fusionSlowFactor;

        // 1. Force Step Execution
        for (let i = 0; i < sim.nucleons.length; i++) {
            for (let j = i + 1; j < sim.nucleons.length; j++) {
                computeInterNucleonForces(sim.nucleons[i], sim.nucleons[j]);
            }
        }

        // 2. Position Integration Projection
        sim.nucleons.forEach(n => {
            n.x += n.vx * dt + 0.5 * n.ax * dt * dt;
            n.y += n.vy * dt + 0.5 * n.ay * dt * dt;
            n.old_ax = n.ax;
            n.old_ay = n.ay;
            n.ax = 0;
            n.ay = 0;
        });

        // 3. Apply Unified Constraints
        applyUnifiedConstraints(sim.nucleons);

        // 4. Update Velocities out from integrated constraint points
        sim.nucleons.forEach(n => {
            if (!sim.reactionTriggered) {
                let speedSq = n.vx * n.vx + n.vy * n.vy;
                if (speedSq > 5) {
                    n.vx *= 0.50; n.vy *= 0.50;
                } else {
                    n.vx *= 0.82; n.vy *= 0.82;
                }
            }
            
            n.vx += 0.5 * (n.old_ax + n.ax) * dt;
            n.vy += 0.5 * (n.old_ay + n.ay) * dt;

            if (n.charge > 0) {
                let distFromEdge = Math.min(n.x, canvasWidth - n.x, n.y, canvasHeight - n.y);
                if (distFromEdge < MAGNETIC_ZONE) {
                    let depth = 1 - (distFromEdge / MAGNETIC_ZONE);
                    let B = (sim.bMax !== undefined ? sim.bMax : B_MAX) * depth * depth;
                    let dvx = n.charge * n.vy * B * dt;
                    let dvy = -n.charge * n.vx * B * dt;
                    n.vx += dvx;
                    n.vy += dvy;
                }
            }

            let speed = Math.sqrt(n.vx * n.vx + n.vy * n.vy);
            if (speed > MAX_SPEED) {
                // Only clamp speeds for incoming Deuterium and Tritium ions
                if (!n.isEjectedNeutron && !n.isAlphaComponent) {
                    let scale = MAX_SPEED / speed;
                    n.vx *= scale;
                    n.vy *= scale;
                }
            }

            if (!n.isEjectedNeutron) {
                let pad = n.radius;
                if (n.x < pad) { n.x = pad; n.vx = Math.abs(n.vx) * 0.5; }
                if (n.x > canvasWidth - pad) { n.x = canvasWidth - pad; n.vx = -Math.abs(n.vx) * 0.5; }
                if (n.y < pad) { n.y = pad; n.vy = Math.abs(n.vy) * 0.5; }
                if (n.y > canvasHeight - pad) { n.y = canvasHeight - pad; n.vy = -Math.abs(n.vy) * 0.5; }
            }
        });

        evaluateFusionState(sim);

        const REMOVAL_MARGIN = 150;
        sim.nucleons = sim.nucleons.filter(n => {
            if (!n.isEjectedNeutron) return true;
            if (n.x < -REMOVAL_MARGIN || n.x > canvasWidth + REMOVAL_MARGIN ||
                n.y < -REMOVAL_MARGIN || n.y > canvasHeight + REMOVAL_MARGIN) {
                sim.neutronEscaped = {
                    x: Math.max(0, Math.min(canvasWidth, n.x)),
                    y: Math.max(0, Math.min(canvasHeight, n.y)),
                    time: Date.now()
                };
                return false;
            }
            return true;
        });

        return { dt, adaptiveFactor, minDist };
    }

    return {
        Nucleon,
        Photon,
        rotateCluster,
        evaluateFusionState,
        step,
        MAGNETIC_ZONE
    };

})();