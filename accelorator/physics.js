const FusionPhysics = (() => {

    const COULOMB_CONSTANT = 1000;
    const MAX_SPEED = 25;

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
            return 0.5 * this.mass * (this.vx * this.vx + this.vy * this.vy);
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

    // Now ONLY computes long-range field effects (Coulomb barrier & initial attraction)
    function computeInterNucleonForces(n1, n2) {
        if (n1.clusterType === n2.clusterType && n1.clusterType !== 'Alpha') return;

        let dx = n2.x - n1.x;
        let dy = n2.y - n1.y;
        let r = Math.sqrt(dx * dx + dy * dy);
        if (r < 0.1) return;

        // 1. Long-range Coulomb Repulsion
        let f_coulomb = (COULOMB_CONSTANT * n1.charge * n2.charge) / (r * r);

        // 2. Long-range Strong Force Attraction
        const MU = 0.15;
        let C_ATTRACT = 1200; 
        if (r < 24.0) C_ATTRACT = 3500; 

        let expTerm = Math.exp(-MU * r);
        let f_strong_attract = C_ATTRACT * (expTerm / r) * (MU + (1 / r));

        // Combined Net macro-force (No Pauli exclusion formula needed anymore!)
        let f_net = f_coulomb - f_strong_attract;

        let fx = (dx / r) * f_net;
        let fy = (dy / r) * f_net;

        n1.ax -= fx / n1.mass;
        n1.ay -= fy / n1.mass;
        n2.ax += fx / n2.mass;
        n2.ay += fy / n2.mass;
    }

    // THE UNIFIED ENGINE CONSTRAINTS: Resolves structural shapes AND overlapping entirely via geometry
    function applyUnifiedConstraints(nucleons) {
        // Run solver passes (More passes = crisper, stiffer physical boundaries)
        const SOLVER_PASSES = 5; 

        for (let pass = 0; pass < SOLVER_PASSES; pass++) {
            for (let i = 0; i < nucleons.length; i++) {
                for (let j = i + 1; j < nucleons.length; j++) {
                    let n1 = nucleons[i];
                    let n2 = nucleons[j];

                    let dx = n2.x - n1.x;
                    let dy = n2.y - n1.y;
                    let dist = Math.sqrt(dx * dx + dy * dy) || 0.1;

                    // Condition A: They belong to the same un-fused parent ion configuration
                    if (n1.clusterType === n2.clusterType && n1.clusterType !== 'Alpha' && n1.clusterType !== 'Free') {
                        let targetDist = (n1.clusterType === 'D') ? 24.0 : 23.5;
                        let difference = targetDist - dist;
                        let percent = (difference / dist) * 0.5;
                        
                        n1.x -= dx * percent; n1.y -= dy * percent;
                        n2.x += dx * percent; n2.y += dy * percent;
                    } 
                    
                    // Condition B: General Solid Boundary Overlapping (Replaces Pauli Exclusion Core)
                    else {
                        let minDist = n1.radius + n2.radius; // 24px hard boundary
                        if (dist < minDist) {
                            let difference = minDist - dist;
                            // Push apart symmetrically based on geometric collision depth
                            let percent = (difference / dist) * 0.5;

                            n1.x -= dx * percent; n1.y -= dy * percent;
                            n2.x += dx * percent; n2.y += dy * percent;
                            
                            // Friction Damping on contact velocities to keep things smooth
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

        let avgX = sim.nucleons.reduce((s, n) => s + n.x, 0) / 5;
        let avgY = sim.nucleons.reduce((s, n) => s + n.y, 0) / 5;

        let maxSpread = Math.max(...sim.nucleons.map(n =>
            Math.sqrt((n.x - avgX) ** 2 + (n.y - avgY) ** 2)
        ));

        // Fusion triggers naturally if macro forces pulled them within capture radius
        if (maxSpread < 42) {
            let neutrons = sim.nucleons.filter(n => n.type === 'neutron');
            neutrons.sort((a, b) =>
                ((b.x - avgX) ** 2 + (b.y - avgY) ** 2) -
                ((a.x - avgX) ** 2 + (a.y - avgY) ** 2)
            );
            let escapeNeutron = neutrons[0];

            sim.nucleons.forEach(n => {
                if (n === escapeNeutron) {
                    let angle = Math.atan2(n.y - avgY, n.x - avgX);
                    n.isEjectedNeutron = true;
                    n.vx = Math.cos(angle) * 16.0;
                    n.vy = Math.sin(angle) * 16.0;
                    n.clusterType = 'Free';
                } else {
                    n.isAlphaComponent = true;
                    n.clusterType = 'Alpha';
                }
            });

            sim.fusionOccurred = true;
            sim.fusionSlowFactor = sim.fusionMinDt;

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

        // 3. Apply Unified Constraints (Instantly repairs any integration errors or overlaps)
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

            let speed = Math.sqrt(n.vx * n.vx + n.vy * n.vy);
            if (speed > MAX_SPEED) {
                let scale = MAX_SPEED / speed;
                n.vx *= scale;
                n.vy *= scale;
            }

            let pad = n.radius;
            if (n.x < pad) { n.x = pad; n.vx *= -0.5; }
            if (n.x > canvasWidth - pad) { n.x = canvasWidth - pad; n.vx *= -0.5; }
            if (n.y < pad) { n.y = pad; n.vy *= -0.5; }
            if (n.y > canvasHeight - pad) { n.y = canvasHeight - pad; n.vy *= -0.5; }
        });

        evaluateFusionState(sim);

        return { dt, adaptiveFactor, minDist };
    }

    return {
        Nucleon,
        Photon,
        rotateCluster,
        evaluateFusionState,
        step
    };

})();