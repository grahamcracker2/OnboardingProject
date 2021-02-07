using Plots
using LinearAlgebra

function vortices()
    v1 = [0,0,0];
    v2 = [0,0,0];
    v3 = [0,0,0];
    v4 = [0,0,0];

    #the current position vectors
    p1 = [0,-0.5,0];
    p2 = [0,0.5,0];
    p3 = [1,0.5,0];
    p4 = [1,-0.5,0];

    backLocation = [0.0];
    frontLocation = [1.0];
    y1 = [-0.5];
    y2 = [0.5];
    y3 = [0.5];
    y4 = [-0.5];

    gamma = [0,0,1];
    dt = 0.01;
    numIterations = 4000;

    #gamma vectors are mulitplied by -1 when they are reversed
    for i = 1:numIterations
        v1 = [0,0,0];
        r12 = p1 - p2; #distance vector from 1 to 2
        r13 = p1 - p3; #distance vector from 1 to 3
        r14 = p1 - p4; #distance vector from 1 to 4
        v1 = v1 + cross(gamma, r12) / (2*pi*norm(r12)^2);
        v1 = v1 + cross(gamma, r13) / (2*pi*norm(r13)^2);
        v1 = v1 + cross(-1 .* gamma, r14) / (2*pi*norm(r14)^2);

        v2 = [0,0,0];
        r21 = p2 - p1;
        r23 = p2 - p3;
        r24 = p2 - p4;
        v2 = v2 + cross(-1 .* gamma, r21) / (2*pi*norm(r21)^2);
        v2 = v2 + cross(gamma, r23) / (2*pi*norm(r23)^2);
        v2 = v2 + cross(-1 .* gamma, r24) / (2*pi*norm(r24)^2);

        v3 = [0,0,0];
        r31 = p3 - p1;
        r32 = p3 - p2;
        r34 = p3 - p4;
        v3 = v3 + cross(-1 .* gamma, r31) / (2*pi*norm(r31)^2);
        v3 = v3 + cross(gamma, r32) / (2*pi*norm(r32)^2);
        v3 = v3 + cross(-1 .* gamma, r34) / (2*pi*norm(r34)^2);

        v4 = [0,0,0];
        r41 = p4 - p1;
        r42 = p4 - p2;
        r43 = p4 - p3;
        v4 = v4 + cross(-1 .* gamma, r41) / (2*pi*norm(r41)^2);
        v4 = v4 + cross(gamma, r42) / (2*pi*norm(r42)^2);
        v4 = v4 + cross(gamma, r43) / (2*pi*norm(r43)^2);

        p1 = p1 + v1.*dt
        p2 = p2 + v2.*dt
        p3 = p3 + v3.*dt
        p4 = p4 + v4.*dt

        append!(backLocation, p1[1]);
        append!(frontLocation, p3[1]);
        append!(y1, p1[2]);
        append!(y2, p2[2]);
        append!(y3, p3[2]);
        append!(y4, p4[2]);
    end
    plot(backLocation, y1, title = "Vortex Ring Progression", legend = :none, linecolor = :red, xlabel = "Lateral Displacement (in)", ylabel = "Vertical Displacement (in)")
    plot!(backLocation, y2, linecolor = :red)
    plot!(frontLocation, y3, linecolor = :black)
    plot!(frontLocation, y4, linecolor = :black)
end

vortices()
