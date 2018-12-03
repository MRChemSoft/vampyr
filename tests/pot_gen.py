def x_tree_generator(x_tree, dir, prec=10e-3, beta=100.0, mid=0.0):
    import vampyr3d as vp
    beta = 100.0
    power = [0, 0, 0]
    power[int(dir)] = 2
    x = vp.GaussFunc(beta, 4.0*beta**2, [mid, mid, mid], power)
    vp.build_grid(x_tree, x)
    vp.project(prec, x_tree, x)


def const_tree_generator(const_tree, prec=10e-3, beta=100.0, mu=0.0, mid=0.0):
    import vampyr3d as vp
    power = [0, 0, 0]
    const = vp.GaussFunc(beta, -6.0*beta - mu**2, [mid, mid, mid], power)
    vp.build_grid(const_tree, const)
    vp.project(prec, const_tree, const)


def v_generator(v_tree, MRA, prec=10e-3, beta=100.0, mid=0.0, mu=0.0):
    import vampyr3d as vp
    from numpy import pi
    alpha = (beta/pi)**(3/2)

    x_tree = vp.FunctionTree(MRA)
    y_tree = vp.FunctionTree(MRA)
    z_tree = vp.FunctionTree(MRA)
    const_tree = vp.FunctionTree(MRA)

    x_tree_generator(x_tree, 0, prec, beta, mid)
    x_tree_generator(y_tree, 1, prec, beta, mid)
    x_tree_generator(z_tree, 2, prec, beta, mid)
    const_tree_generator(const_tree, prec, beta, mu, mid)
    vec = []
    vec.append(tuple([1.0, x_tree]))
    vec.append(tuple([1.0, y_tree]))
    vec.append(tuple([1.0, z_tree]))
    vec.append(tuple([1.0, const_tree]))
    vp.build_grid(v_tree, vec)
    vp.add(prec, v_tree, vec)

    alpha = -alpha/(4.0*pi)
    v_tree.rescale(alpha)
