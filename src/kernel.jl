module D2Q9
    const celerities = transpose(Int32[
        0 0;
        1 0; 0 1; -1 0; 0 -1;
        1 1; -1 1; -1 -1; 1 -1
    ])
    const weights = vcat(4./9., [1./9. for i=1:4], [1./36. for i=1:4])
    const speed_of_sound = 1./sqrt(3.)
end

module D3Q19
    const celerities = transpose(Int32[
        0 0 0;
        1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1;
        1 1 0; 1 0 1; -1 1 0; -1 0 1; -1 -1 0; -1 0 -1; 1 -1 0; 1 0 -1;
        0 1 1; 0 -1 1; 0 -1 -1; 0 1 -1;
    ])
    const weights = vcat(1./3., [1./18. for i=1:6], [7./18. for i=7:18])
    const speed_of_sound = (
        sum([weights[i] * vecnorm(celerities) for i=1:length(weights)])
    )
end
