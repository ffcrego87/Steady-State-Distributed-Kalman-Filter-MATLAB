function noise = gen_noise(dim,bound)
    noise=bound*(2*(rand(dim,1)>0.5)-1);
end