function kernel_eps = get_kernel_width(N)

    assert(N>0);    
    kernel_eps = log(N)^(3/8)/N^(1/4);

end