import numpy as np

####################################################################################################


def RemoteNodeForce(x1, x2, x3, x4, bp, b, a, mu, nu):
    # this calculates the forces between dislocation nodes analytically
    # inputs: endpoints of first dislocation segment starting at x1 ending at x2 with burgers vector bp
    #        endpoints of second dislocation segment starting at x3 ending at x4 with burgers vector b
    #        core paramter a
    #        shear modulus mu
    #        poisson ration nu
    #
    # outputs: f1,f2,f3,f4 is the force on nodes located at x1, x2, x3, x4 respectively
    if x1.ndim == 1:
        x1 = np.expand_dims(x1, axis=0)
        x2 = np.expand_dims(x2, axis=0)
        x3 = np.expand_dims(x3, axis=0)
        x4 = np.expand_dims(x4, axis=0)
        bp = np.expand_dims(bp, axis=0)
        b = np.expand_dims(b, axis=0)

    f1 = np.zeros_like(x1)
    f2 = np.zeros_like(x2)
    f4 = np.zeros_like(x3)
    f3 = np.zeros_like(x4)
    eps = 1e-06

    Diff = x4 - x3
    oneoverL = 1.0 / np.sqrt(np.sum(np.multiply(Diff, Diff), axis=1))
    t = np.multiply(Diff, np.array([oneoverL, oneoverL, oneoverL]).T)
    Diff = x2 - x1
    oneoverLp = 1.0 / np.sqrt(np.sum(np.multiply(Diff, Diff), axis=1))
    tp = np.multiply(Diff, np.array([oneoverLp, oneoverLp, oneoverLp]).T)

    c = np.sum(np.multiply(t, tp), axis=1)
    c2 = np.multiply(c, c)
    onemc2 = 1 - c2

    # Find the indices of the segments that are parallel (need special treatment)
    spindex = onemc2 < eps
    x1sp = np.copy(x1[spindex, :])
    x2sp = np.copy(x2[spindex, :])
    x3sp = np.copy(x3[spindex, :])
    x4sp = np.copy(x4[spindex, :])
    bsp = np.copy(b[spindex, :])
    bpsp = np.copy(bp[spindex, :])

    nonpar = np.logical_not(spindex)
    if np.sum(nonpar) > 0:
        txtp = np.cross(t[nonpar, :], tp[nonpar, :])
        onemc2inv = 1.0 / onemc2[nonpar]
        R = np.concatenate(
            (x3[nonpar, :] - x1[nonpar, :], x4[nonpar, :] - x2[nonpar, :]), axis=1
        )
        d = np.multiply(np.sum(np.multiply(R[:, 0:3], txtp), 1), onemc2inv)
        temp1 = np.array(
            [np.sum(np.multiply(R[:, 0:3], t), 1), np.sum(np.multiply(R[:, 3:6], t), 1)]
        ).T
        temp2 = np.array(
            [
                np.sum(np.multiply(R[:, 0:3], tp), 1),
                np.sum(np.multiply(R[:, 3:6], tp), 1),
            ]
        ).T
        y = np.multiply(
            (temp1 - np.multiply(np.array([c, c]).transpose(), temp2)),
            np.array([onemc2inv, onemc2inv]).T,
        )
        z = np.multiply(
            (temp2 - np.multiply(np.array([c, c]).transpose(), temp1)),
            np.array([onemc2inv, onemc2inv]).T,
        )
        yin = np.array([y[:, 0], y[:, 0], y[:, 1], y[:, 1]]).T
        zin = np.array([z[:, 0], z[:, 1], z[:, 0], z[:, 1]]).T

        ####################################################################################################################
        #  this section calculates the formulae from the integral expressions
        a2 = a * a
        a2_d2 = a2 + np.multiply(np.multiply(d, d), onemc2)
        y2 = np.multiply(yin, yin)
        z2 = np.multiply(zin, zin)
        Ra = np.sqrt(
            np.array([a2_d2, a2_d2, a2_d2, a2_d2]).T
            + y2
            + z2
            + np.multiply(np.multiply(2.0 * yin, zin), np.array([c, c, c, c]).T)
        )
        Rainv = 1.0 / Ra

        Ra_Rdot_tp = Ra + zin + np.multiply(yin, np.array([c, c, c, c]).T)
        Ra_Rdot_t = Ra + yin + np.multiply(zin, np.array([c, c, c, c]).T)

        log_Ra_Rdot_tp = np.log(Ra_Rdot_tp)
        ylog_Ra_Rdot_tp = np.multiply(yin, log_Ra_Rdot_tp)

        log_Ra_Rdot_t = np.log(Ra_Rdot_t)
        zlog_Ra_Rdot_t = np.multiply(zin, log_Ra_Rdot_t)

        Ra2_R_tpinv = Rainv / Ra_Rdot_tp
        yRa2_R_tpinv = np.multiply(yin, Ra2_R_tpinv)
        y2Ra2_R_tpinv = np.multiply(yin, yRa2_R_tpinv)

        Ra2_R_tinv = Rainv / Ra_Rdot_t
        zRa2_R_tinv = np.multiply(zin, Ra2_R_tinv)
        z2Ra2_R_tinv = np.multiply(zin, zRa2_R_tinv)

        denom = 1.0 / np.sqrt(np.multiply(onemc2, a2_d2))
        cdenom = np.multiply((1 + c), denom)

        f_003 = np.multiply(
            -2.0 * np.array([denom, denom, denom, denom]).T,
            np.arctan(
                np.multiply(
                    (Ra + yin + zin), np.array([cdenom, cdenom, cdenom, cdenom]).T
                )
            ),
        )
        adf_003 = np.multiply(np.array([a2_d2, a2_d2, a2_d2, a2_d2]).T, f_003)
        commonf223 = np.multiply(
            (np.multiply(np.array([c, c, c, c]).T, Ra) - adf_003),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )

        f_103 = np.multiply(
            (np.multiply(np.array([c, c, c, c]).T, log_Ra_Rdot_t) - log_Ra_Rdot_tp),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_013 = np.multiply(
            (np.multiply(np.array([c, c, c, c]).T, log_Ra_Rdot_tp) - log_Ra_Rdot_t),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_113 = np.multiply(
            (np.multiply(np.array([c, c, c, c]).T, adf_003) - Ra),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_203 = zlog_Ra_Rdot_t + commonf223
        f_023 = ylog_Ra_Rdot_tp + commonf223

        commonf225 = f_003 - np.multiply(np.array([c, c, c, c]).T, Rainv)
        commonf025 = np.multiply(np.array([c, c, c, c]).T, yRa2_R_tpinv) - Rainv
        ycommonf025 = np.multiply(yin, commonf025)
        commonf205 = np.multiply(np.array([c, c, c, c]).T, zRa2_R_tinv) - Rainv
        zcommonf205 = np.multiply(zin, commonf205)
        commonf305 = (
            log_Ra_Rdot_t
            - np.multiply((yin - np.multiply(np.array([c, c, c, c]).T, zin)), Rainv)
            - np.multiply(np.array([c2, c2, c2, c2]).T, z2Ra2_R_tinv)
        )
        zcommonf305 = np.multiply(zin, commonf305)
        commonf035 = (
            log_Ra_Rdot_tp
            - np.multiply((zin - np.multiply(np.array([c, c, c, c]).T, yin)), Rainv)
            - np.multiply(np.array([c2, c2, c2, c2]).T, y2Ra2_R_tpinv)
        )
        tf_113 = 2.0 * f_113

        f_005 = (f_003 - yRa2_R_tpinv - zRa2_R_tinv) / np.array(
            [a2_d2, a2_d2, a2_d2, a2_d2]
        ).T
        f_105 = np.multiply(
            (Ra2_R_tpinv - np.multiply(np.array([c, c, c, c]).T, Ra2_R_tinv)),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_015 = np.multiply(
            (Ra2_R_tinv - np.multiply(np.array([c, c, c, c]).T, Ra2_R_tpinv)),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_115 = np.multiply(
            (
                Rainv
                - np.multiply(
                    np.array([c, c, c, c]).T, (yRa2_R_tpinv + zRa2_R_tinv + f_003)
                )
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_205 = np.multiply(
            (
                yRa2_R_tpinv
                + np.multiply(np.array([c2, c2, c2, c2]).T, zRa2_R_tinv)
                + commonf225
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_025 = np.multiply(
            (
                zRa2_R_tinv
                + np.multiply(np.array([c2, c2, c2, c2]).T, yRa2_R_tpinv)
                + commonf225
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_215 = np.multiply(
            (
                f_013
                - ycommonf025
                + np.multiply(np.array([c, c, c, c]).T, (zcommonf205 - f_103))
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_125 = np.multiply(
            (
                f_103
                - zcommonf205
                + np.multiply(np.array([c, c, c, c]).T, (ycommonf025 - f_013))
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_225 = np.multiply(
            (
                f_203
                - zcommonf305
                + np.multiply(
                    np.array([c, c, c, c]).T, (np.multiply(y2, commonf025) - tf_113)
                )
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_305 = np.multiply(
            (
                y2Ra2_R_tpinv
                + np.multiply(np.array([c, c, c, c]).T, commonf305)
                + 2.0 * f_103
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_035 = np.multiply(
            (
                z2Ra2_R_tinv
                + np.multiply(np.array([c, c, c, c]).T, commonf035)
                + 2.0 * f_013
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_315 = np.multiply(
            (
                tf_113
                - np.multiply(y2, commonf025)
                + np.multiply(np.array([c, c, c, c]).T, (zcommonf305 - f_203))
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        f_135 = np.multiply(
            (
                tf_113
                - np.multiply(z2, commonf205)
                + np.multiply(
                    np.array([c, c, c, c]).T, (np.multiply(yin, commonf035) - f_023)
                )
            ),
            np.array([onemc2inv, onemc2inv, onemc2inv, onemc2inv]).T,
        )
        mf = np.array([1, -1, -1, 1])
        Fintegrals = np.array(
            [
                f_003 @ mf,
                f_103 @ mf,
                f_013 @ mf,
                f_113 @ mf,
                f_203 @ mf,
                f_023 @ mf,
                f_005 @ mf,
                f_105 @ mf,
                f_015 @ mf,
                f_115 @ mf,
                f_205 @ mf,
                f_025 @ mf,
                f_215 @ mf,
                f_125 @ mf,
                f_225 @ mf,
                f_305 @ mf,
                f_035 @ mf,
                f_315 @ mf,
                f_135 @ mf,
            ]
        ).T

        #################################################################################################################
        # this section calculates the dot products and cross prodcucts for the coefficients
        m4p = 0.25 * mu / np.pi
        m4pd = np.multiply(m4p, d)
        m8p = 0.5 * m4p
        m8pd = np.multiply(m8p, d)
        m4pn = m4p / (1 - nu)
        m4pnd = np.multiply(m4pn, d)
        m4pnd2 = np.multiply(m4pnd, d)
        m4pnd3 = np.multiply(m4pnd2, d)
        a2m4pnd = np.multiply(a2, m4pnd)
        a2m8pd = np.multiply(a2, m8pd)
        a2m4pn = a2 * m4pn
        a2m8p = a2 * m8p

        tpxt = -txtp
        txbp = np.cross(t[nonpar, :], bp[nonpar, :])
        tpxb = np.cross(tp[nonpar, :], b[nonpar, :])
        bxt = np.cross(b[nonpar, :], t[nonpar, :])
        bpxtp = np.cross(bp[nonpar, :], tp[nonpar, :])

        tdb = np.sum(np.multiply(t[nonpar, :], b[nonpar, :]), 1)
        tdbp = np.sum(np.multiply(t[nonpar, :], bp[nonpar, :]), 1)
        tpdb = np.sum(np.multiply(tp[nonpar, :], b[nonpar, :]), 1)
        tpdbp = np.sum(np.multiply(tp[nonpar, :], bp[nonpar, :]), 1)
        txtpdb = np.sum(np.multiply(txtp[nonpar, :], b[nonpar, :]), 1)
        tpxtdbp = np.sum(np.multiply(tpxt[nonpar, :], bp[nonpar, :]), 1)
        txbpdtp = tpxtdbp
        tpxbdt = txtpdb

        bpxtpdb = np.sum(np.multiply(bpxtp, b), 1)
        bxtdbp = np.sum(np.multiply(bxt, bp), 1)
        txbpdb = bxtdbp
        tpxbdbp = bpxtpdb
        txtpxt = tp - np.multiply(np.array([c, c, c]).T, t)
        tpxtxtp = t - np.multiply(np.array([c, c, c]).T, tp)
        txtpxbp = np.multiply(np.array([tdbp, tdbp, tdbp]).T, tp) - np.multiply(
            np.array([tpdbp, tpdbp, tpdbp]).T, t
        )
        tpxtxb = np.multiply(np.array([tpdb, tpdb, tpdb]).T, t) - np.multiply(
            np.array([tdb, tdb, tdb]).T, tp
        )
        txbpxt = bp - np.multiply(np.array([tdbp, tdbp, tdbp]).T, t)
        tpxbxtp = b - np.multiply(np.array([tpdb, tpdb, tpdb]).T, tp)
        bpxtpxt = np.multiply(np.array([tdbp, tdbp, tdbp]).T, tp) - np.multiply(
            np.array([c, c, c]).T, bp
        )
        bxtxtp = np.multiply(np.array([tpdb, tpdb, tpdb]).T, t) - np.multiply(
            np.array([c, c, c]).T, b
        )
        txtpxbpxt = np.multiply(np.array([tdbp, tdbp, tdbp]).T, tpxt)
        tpxtxbxtp = np.multiply(np.array([tpdb, tpdb, tpdb]).T, txtp)
        txtpxbpdtp = tdbp - np.multiply(tpdbp, c)
        tpxtxbdt = tpdb - np.multiply(tdb, c)
        txtpxbpdb = np.multiply(tdbp, tpdb) - np.multiply(tpdbp, tdb)
        tpxtxbdbp = txtpxbpdb

        ##################################################################################################################
        # this section calculates the coefficients for two of the forces
        temp1 = np.multiply(tdbp, tpdb) + txtpxbpdb
        I00a = np.multiply(np.array([temp1, temp1, temp1]).T, tpxt)
        I00b = np.multiply(bxt, np.array([txtpxbpdtp, txtpxbpdtp, txtpxbpdtp]).T)
        temp1 = np.multiply(m4pnd, txtpdb)
        temp2 = np.multiply(m4pnd, bpxtpdb)
        I_003 = (
            np.multiply(np.array([m4pd, m4pd, m4pd]).T, I00a)
            - np.multiply(np.array([m4pnd, m4pnd, m4pnd]).T, I00b)
            + np.multiply(np.array([temp1, temp1, temp1]).T, bpxtpxt)
            + np.multiply(np.array([temp2, temp2, temp2]).T, txtpxt)
        )
        temp1 = np.multiply(np.multiply(m4pnd3, txtpxbpdtp), txtpdb)
        I_005 = (
            np.multiply(np.array([a2m8pd, a2m8pd, a2m8pd]).T, I00a)
            - np.multiply(np.array([a2m4pnd, a2m4pnd, a2m4pnd]).T, I00b)
            - np.multiply(np.array([temp1, temp1, temp1]).T, txtpxt)
        )
        I10a = np.multiply(txbpxt, np.array([tpdb, tpdb, tpdb]).T) - np.multiply(
            txtp, np.array([txbpdb, txbpdb, txbpdb]).T
        )
        I10b = np.multiply(bxt, np.array([txbpdtp, txbpdtp, txbpdtp]).T)
        temp1 = np.multiply(m4pn, tdb)
        I_103 = (
            np.multiply(np.array([temp1, temp1, temp1]).T, bpxtpxt)
            + np.multiply(m4p, I10a)
            - np.multiply(m4pn, I10b)
        )
        temp1 = np.multiply(
            m4pnd2, (np.multiply(txbpdtp, txtpdb) + np.multiply(txtpxbpdtp, tdb))
        )
        I_105 = (
            np.multiply(a2m8p, I10a)
            - np.multiply(a2m4pn, I10b)
            - np.multiply(np.array([temp1, temp1, temp1]).T, txtpxt)
        )
        I01a = np.multiply(txtp, np.array([bpxtpdb, bpxtpdb, bpxtpdb]).T) - np.multiply(
            bpxtpxt, np.array([tpdb, tpdb, tpdb]).T
        )
        temp1 = np.multiply(m4pn, tpdb)
        temp2 = np.multiply(m4pn, bpxtpdb)
        I_013 = (
            np.multiply(m4p, I01a)
            + np.multiply(np.array([temp1, temp1, temp1]).T, bpxtpxt)
            - np.multiply(np.array([temp2, temp2, temp2]).T, txtp)
        )
        temp1 = np.multiply(np.multiply(m4pnd2, txtpxbpdtp), tpdb)
        temp2 = np.multiply(np.multiply(m4pnd2, txtpxbpdtp), txtpdb)
        I_015 = (
            np.multiply(a2m8p, I01a)
            - np.multiply(np.array([temp1, temp1, temp1]).T, txtpxt)
            + np.multiply(np.array([temp2, temp2, temp2]).T, txtp)
        )
        temp1 = np.multiply(np.multiply(m4pnd, txbpdtp), tdb)
        I_205 = np.multiply(-np.array([temp1, temp1, temp1]).T, txtpxt)
        temp1 = np.multiply(np.multiply(m4pnd, txtpxbpdtp), tpdb)
        I_025 = np.multiply(np.array([temp1, temp1, temp1]).T, txtp)
        temp1 = np.multiply(
            m4pnd, (np.multiply(txtpxbpdtp, tdb) + np.multiply(txbpdtp, txtpdb))
        )
        temp2 = np.multiply(np.multiply(m4pnd, txbpdtp), tpdb)
        I_115 = np.multiply(np.array([temp1, temp1, temp1]).T, txtp) - np.multiply(
            np.array([temp2, temp2, temp2]).T, txtpxt
        )
        temp1 = np.multiply(np.multiply(m4pn, txbpdtp), tdb)
        I_215 = np.multiply(np.array([temp1, temp1, temp1]).T, txtp)
        temp1 = np.multiply(np.multiply(m4pn, txbpdtp), tpdb)
        I_125 = np.multiply(np.array([temp1, temp1, temp1]).T, txtp)

        ####################################################################################################################
        # this section calculates the first two forces
        Fint_003 = Fintegrals[:, 2 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 1 - 1])
        Fint_103 = Fintegrals[:, 5 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 2 - 1])
        Fint_013 = Fintegrals[:, 4 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 3 - 1])
        Fint_005 = Fintegrals[:, 8 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 7 - 1])
        Fint_105 = Fintegrals[:, 11 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 8 - 1]
        )
        Fint_015 = Fintegrals[:, 10 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 9 - 1]
        )
        Fint_115 = Fintegrals[:, 13 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 10 - 1]
        )
        Fint_205 = Fintegrals[:, 16 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 11 - 1]
        )
        Fint_025 = Fintegrals[:, 14 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 12 - 1]
        )
        Fint_215 = Fintegrals[:, 18 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 13 - 1]
        )
        Fint_125 = Fintegrals[:, 15 - 1] - np.multiply(
            y[:, 1 - 1], Fintegrals[:, 14 - 1]
        )
        f4 = (
            np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T)
            + np.multiply(I_103, np.array([Fint_103, Fint_103, Fint_103]).T)
            + np.multiply(I_013, np.array([Fint_013, Fint_013, Fint_013]).T)
        )
        f4 = (
            f4
            + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
            + np.multiply(I_105, np.array([Fint_105, Fint_105, Fint_105]).T)
            + np.multiply(I_015, np.array([Fint_015, Fint_015, Fint_015]).T)
        )
        f4 = (
            f4
            + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
            + np.multiply(I_205, np.array([Fint_205, Fint_205, Fint_205]).T)
            + np.multiply(I_025, np.array([Fint_025, Fint_025, Fint_025]).T)
        )
        f4 = (
            f4
            + np.multiply(I_215, np.array([Fint_215, Fint_215, Fint_215]).T)
            + np.multiply(I_125, np.array([Fint_125, Fint_125, Fint_125]).T)
        )
        f4 = np.multiply(f4, np.array([oneoverL, oneoverL, oneoverL]).T)

        Fint_003 = np.multiply(y[:, 2 - 1], Fintegrals[:, 1 - 1]) - Fintegrals[:, 2 - 1]
        Fint_103 = np.multiply(y[:, 2 - 1], Fintegrals[:, 2 - 1]) - Fintegrals[:, 5 - 1]
        Fint_013 = np.multiply(y[:, 2 - 1], Fintegrals[:, 3 - 1]) - Fintegrals[:, 4 - 1]
        Fint_005 = np.multiply(y[:, 2 - 1], Fintegrals[:, 7 - 1]) - Fintegrals[:, 8 - 1]
        Fint_105 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 8 - 1]) - Fintegrals[:, 11 - 1]
        )
        Fint_015 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 9 - 1]) - Fintegrals[:, 10 - 1]
        )
        Fint_115 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 10 - 1]) - Fintegrals[:, 13 - 1]
        )
        Fint_205 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 11 - 1]) - Fintegrals[:, 16 - 1]
        )
        Fint_025 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 12 - 1]) - Fintegrals[:, 14 - 1]
        )
        Fint_215 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 13 - 1]) - Fintegrals[:, 18 - 1]
        )
        Fint_125 = (
            np.multiply(y[:, 2 - 1], Fintegrals[:, 14 - 1]) - Fintegrals[:, 15 - 1]
        )
        f3 = (
            np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T)
            + np.multiply(I_103, np.array([Fint_103, Fint_103, Fint_103]).T)
            + np.multiply(I_013, np.array([Fint_013, Fint_013, Fint_013]).T)
        )
        f3 = (
            f3
            + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
            + np.multiply(I_105, np.array([Fint_105, Fint_105, Fint_105]).T)
            + np.multiply(I_015, np.array([Fint_015, Fint_015, Fint_015]).T)
        )
        f3 = (
            f3
            + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
            + np.multiply(I_205, np.array([Fint_205, Fint_205, Fint_205]).T)
            + np.multiply(I_025, np.array([Fint_025, Fint_025, Fint_025]).T)
        )
        f3 = (
            f3
            + np.multiply(I_215, np.array([Fint_215, Fint_215, Fint_215]).T)
            + np.multiply(I_125, np.array([Fint_125, Fint_125, Fint_125]).T)
        )
        f3 = np.multiply(f3, np.array([oneoverL, oneoverL, oneoverL]).T)

        #######################################################################################################################
        # this section calculates the coefficients for the second two forces
        temp1 = np.multiply(tpdb, tdbp) + tpxtxbdbp
        I00a = np.multiply(np.array([temp1, temp1, temp1]).T, txtp)
        I00b = np.multiply(bpxtp, np.array([tpxtxbdt, tpxtxbdt, tpxtxbdt]).T)
        temp1 = np.multiply(m4pnd, tpxtdbp)
        temp2 = np.multiply(m4pnd, bxtdbp)
        I_003 = (
            np.multiply(np.array([m4pd, m4pd, m4pd]).T, I00a)
            - np.multiply(np.array([m4pnd, m4pnd, m4pnd]).T, I00b)
            + np.multiply(np.array([temp1, temp1, temp1]).T, bxtxtp)
            + np.multiply(np.array([temp2, temp2, temp2]).T, tpxtxtp)
        )
        temp1 = np.multiply(np.multiply(m4pnd3, tpxtxbdt), tpxtdbp)
        I_005 = (
            np.multiply(np.array([a2m8pd, a2m8pd, a2m8pd]).T, I00a)
            - np.multiply(np.array([a2m4pnd, a2m4pnd, a2m4pnd]).T, I00b)
            - np.multiply(np.array([temp1, temp1, temp1]).T, tpxtxtp)
        )
        I01a = np.multiply(tpxt, np.array([tpxbdbp, tpxbdbp, tpxbdbp]).T) - np.multiply(
            tpxbxtp, np.array([tdbp, tdbp, tdbp]).T
        )
        I01b = np.multiply(-bpxtp, np.array([tpxbdt, tpxbdt, tpxbdt]).T)
        temp1 = np.multiply(m4pn, tpdbp)
        I_013 = (
            np.multiply(-np.array([temp1, temp1, temp1]).T, bxtxtp)
            + np.multiply(m4p, I01a)
            - np.multiply(m4pn, I01b)
        )
        temp1 = np.multiply(
            m4pnd2, (np.multiply(tpxbdt, tpxtdbp) + np.multiply(tpxtxbdt, tpdbp))
        )
        I_015 = (
            np.multiply(a2m8p, I01a)
            - np.multiply(a2m4pn, I01b)
            + np.multiply(np.array([temp1, temp1, temp1]).T, tpxtxtp)
        )
        I10a = np.multiply(bxtxtp, np.array([tdbp, tdbp, tdbp]).T) - np.multiply(
            tpxt, np.array([bxtdbp, bxtdbp, bxtdbp]).T
        )
        temp1 = np.multiply(m4pn, tdbp)
        temp2 = np.multiply(m4pn, bxtdbp)
        I_103 = (
            np.multiply(m4p, I10a)
            - np.multiply(np.array([temp1, temp1, temp1]).T, bxtxtp)
            + np.multiply(np.array([temp2, temp2, temp2]).T, tpxt)
        )
        temp1 = np.multiply(np.multiply(m4pnd2, tpxtxbdt), tdbp)
        temp2 = np.multiply(np.multiply(m4pnd2, tpxtxbdt), tpxtdbp)
        I_105 = (
            np.multiply(a2m8p, I10a)
            + np.multiply(np.array([temp1, temp1, temp1]).T, tpxtxtp)
            - np.multiply(np.array([temp2, temp2, temp2]).T, tpxt)
        )
        temp1 = np.multiply(np.multiply(m4pnd, tpxbdt), tpdbp)
        I_025 = np.multiply(-np.array([temp1, temp1, temp1]).T, tpxtxtp)
        temp1 = np.multiply(np.multiply(m4pnd, tpxtxbdt), tdbp)
        I_205 = np.multiply(np.array([temp1, temp1, temp1]).T, tpxt)
        temp1 = np.multiply(
            m4pnd, (np.multiply(tpxtxbdt, tpdbp) + np.multiply(tpxbdt, tpxtdbp))
        )
        temp2 = np.multiply(np.multiply(m4pnd, tpxbdt), tdbp)
        I_115 = np.multiply(np.array([temp1, temp1, temp1]).T, tpxt) - np.multiply(
            np.array([temp2, temp2, temp2]).T, tpxtxtp
        )
        temp1 = np.multiply(np.multiply(m4pn, tpxbdt), tpdbp)
        I_125 = np.multiply(-np.array([temp1, temp1, temp1]).T, tpxt)
        temp1 = np.multiply(np.multiply(m4pn, tpxbdt), tdbp)
        I_215 = np.multiply(-np.array([temp1, temp1, temp1]).T, tpxt)
        ########################################################################################################################
        # this section calculates the second two forces
        Fint_003 = Fintegrals[:, 3 - 1] - np.multiply(z[:, 2 - 1], Fintegrals[:, 1 - 1])
        Fint_103 = Fintegrals[:, 4 - 1] - np.multiply(z[:, 2 - 1], Fintegrals[:, 2 - 1])
        Fint_013 = Fintegrals[:, 6 - 1] - np.multiply(z[:, 2 - 1], Fintegrals[:, 3 - 1])
        Fint_005 = Fintegrals[:, 9 - 1] - np.multiply(z[:, 2 - 1], Fintegrals[:, 7 - 1])
        Fint_105 = Fintegrals[:, 10 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 8 - 1]
        )
        Fint_015 = Fintegrals[:, 12 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 9 - 1]
        )
        Fint_115 = Fintegrals[:, 14 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 10 - 1]
        )
        Fint_205 = Fintegrals[:, 13 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 11 - 1]
        )
        Fint_025 = Fintegrals[:, 17 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 12 - 1]
        )
        Fint_215 = Fintegrals[:, 15 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 13 - 1]
        )
        Fint_125 = Fintegrals[:, 19 - 1] - np.multiply(
            z[:, 2 - 1], Fintegrals[:, 14 - 1]
        )
        f1 = (
            np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T)
            + np.multiply(I_103, np.array([Fint_103, Fint_103, Fint_103]).T)
            + np.multiply(I_013, np.array([Fint_013, Fint_013, Fint_013]).T)
        )
        f1 = (
            f1
            + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
            + np.multiply(I_105, np.array([Fint_105, Fint_105, Fint_105]).T)
            + np.multiply(I_015, np.array([Fint_015, Fint_015, Fint_015]).T)
        )
        f1 = (
            f1
            + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
            + np.multiply(I_205, np.array([Fint_205, Fint_205, Fint_205]).T)
            + np.multiply(I_025, np.array([Fint_025, Fint_025, Fint_025]).T)
        )
        f1 = (
            f1
            + np.multiply(I_215, np.array([Fint_215, Fint_215, Fint_215]).T)
            + np.multiply(I_125, np.array([Fint_125, Fint_125, Fint_125]).T)
        )
        f1 = np.multiply(f1, np.array([oneoverLp, oneoverLp, oneoverLp]).T)
        Fint_003 = np.multiply(z[:, 1 - 1], Fintegrals[:, 1 - 1]) - Fintegrals[:, 3 - 1]
        Fint_103 = np.multiply(z[:, 1 - 1], Fintegrals[:, 2 - 1]) - Fintegrals[:, 4 - 1]
        Fint_013 = np.multiply(z[:, 1 - 1], Fintegrals[:, 3 - 1]) - Fintegrals[:, 6 - 1]
        Fint_005 = np.multiply(z[:, 1 - 1], Fintegrals[:, 7 - 1]) - Fintegrals[:, 9 - 1]
        Fint_105 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 8 - 1]) - Fintegrals[:, 10 - 1]
        )
        Fint_015 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 9 - 1]) - Fintegrals[:, 12 - 1]
        )
        Fint_115 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 10 - 1]) - Fintegrals[:, 14 - 1]
        )
        Fint_205 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 11 - 1]) - Fintegrals[:, 13 - 1]
        )
        Fint_025 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 12 - 1]) - Fintegrals[:, 17 - 1]
        )
        Fint_215 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 13 - 1]) - Fintegrals[:, 15 - 1]
        )
        Fint_125 = (
            np.multiply(z[:, 1 - 1], Fintegrals[:, 14 - 1]) - Fintegrals[:, 19 - 1]
        )
        f2 = (
            np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T)
            + np.multiply(I_103, np.array([Fint_103, Fint_103, Fint_103]).T)
            + np.multiply(I_013, np.array([Fint_013, Fint_013, Fint_013]).T)
        )
        f2 = (
            f2
            + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
            + np.multiply(I_105, np.array([Fint_105, Fint_105, Fint_105]).T)
            + np.multiply(I_015, np.array([Fint_015, Fint_015, Fint_015]).T)
        )
        f2 = (
            f2
            + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
            + np.multiply(I_205, np.array([Fint_205, Fint_205, Fint_205]).T)
            + np.multiply(I_025, np.array([Fint_025, Fint_025, Fint_025]).T)
        )
        f2 = (
            f2
            + np.multiply(I_215, np.array([Fint_215, Fint_215, Fint_215]).T)
            + np.multiply(I_125, np.array([Fint_125, Fint_125, Fint_125]).T)
        )
        f2 = np.multiply(f2, np.array([oneoverLp, oneoverLp, oneoverLp]).T)

    if np.sum(spindex) > 0:
        # this is the parallel case the two lines are parallel use a special lower dimensional function
        f1sp, f2sp, f3sp, f4sp = SpecialRemoteNodeForce(
            x1sp, x2sp, x3sp, x4sp, bpsp, bsp, a, mu, nu, eps
        )
        f1[spindex, :] = f1sp
        f2[spindex, :] = f2sp
        f3[spindex, :] = f3sp
        f4[spindex, :] = f4sp

    return f1, f2, f3, f4


####################################################################################################


def SpecialRemoteNodeForce(
    x1=None,
    x2=None,
    x3=None,
    x4=None,
    bp=None,
    b=None,
    a=None,
    mu=None,
    nu=None,
    ecrit=None,
):
    # this calculates the forces between dislocation nodes analytically
    # this is a special subroutine used for dislocation segments that are too close to parallel to be
    # calculated by the regular expression for forces
    # inputs: endpoints of first dislocation segment starting at x1 ending at x2 with burgers vector bp
    #        endpoints of second dislocation segment starting at x3 ending at x4 with burgers vector b
    #        core paramter a
    #        shear modulus mu
    #        poisson ration nu
    #
    # outputs: f1,f2,f3,f4 is the force on nodes located at x1, x2, x3, x4 respectively

    cotanthetac = np.sqrt((1 - ecrit * (1.01)) / (ecrit * 1.01))
    eps = 1e-16
    Diff = x4 - x3
    oneoverL = 1.0 / np.sqrt(np.sum(np.multiply(Diff, Diff), axis=1))
    t = np.multiply(Diff, np.array([oneoverL, oneoverL, oneoverL]).T)

    Diff = x2 - x1
    oneoverLp = 1.0 / np.sqrt(np.sum(np.multiply(Diff, Diff), axis=1))
    tp = np.multiply(Diff, np.array([oneoverLp, oneoverLp, oneoverLp]).T)
    c = np.sum(np.multiply(t, tp), axis=1)

    flip = c < 0
    x1_old = np.copy(x1)
    x2_old = np.copy(x2)
    x1[flip, :] = x2_old[flip, :]
    x2[flip, :] = x1_old[flip, :]
    tp[flip, :] = -tp[flip, :]
    bp[flip, :] = -bp[flip, :]

    temp = np.sum(np.multiply((x2 - x1), t), axis=1)
    x2mod = x1 + np.multiply(np.array([temp, temp, temp]).T, t)
    diff = x2 - x2mod
    magdiff = np.sqrt(np.sum(np.multiply(diff, diff), axis=1))
    temp = np.multiply(
        np.multiply((0.5 * cotanthetac), np.array([magdiff, magdiff, magdiff])).T, t
    )
    x1mod = x1 + np.multiply(0.5, diff) + temp
    x2mod = x2mod + np.multiply(0.5, diff) - temp
    R = x3 - x1mod
    Rdt = np.sum(np.multiply(R, t), axis=1)
    nd = R - np.multiply(np.array([Rdt, Rdt, Rdt]).T, t)
    d2 = np.sum(np.multiply(nd, nd), axis=1)

    r4 = np.sum(np.multiply(x4, t), axis=1)
    r3 = np.sum(np.multiply(x3, t), axis=1)
    s2 = np.sum(np.multiply(x2mod, t), axis=1)
    s1 = np.sum(np.multiply(x1mod, t), axis=1)

    y = np.array([r3, r3, r4, r4]).T
    z = np.array([-s1, -s2, -s1, -s2]).T

    a2 = a * a
    a2_d2 = a2 + d2
    temp = 1.0 / a2_d2
    a2d2inv = np.array([temp, temp, temp, temp]).T
    ypz = y + z
    ymz = y - z
    Ra = np.sqrt(np.array([a2_d2, a2_d2, a2_d2, a2_d2]).T + np.multiply(ypz, ypz))
    Rainv = 1.0 / Ra
    Log_Ra_ypz = np.log(Ra + ypz)

    f_003 = np.multiply(Ra, a2d2inv)
    f_103 = np.multiply(-0.5, (Log_Ra_ypz - np.multiply(np.multiply(ymz, Ra), a2d2inv)))
    f_013 = np.multiply(-0.5, (Log_Ra_ypz + np.multiply(np.multiply(ymz, Ra), a2d2inv)))
    f_113 = -Log_Ra_ypz
    f_213 = np.multiply(z, Log_Ra_ypz) - Ra
    f_123 = np.multiply(y, Log_Ra_ypz) - Ra

    f_005 = np.multiply(a2d2inv, (np.multiply(2.0 * a2d2inv, Ra) - Rainv))
    f_105 = np.multiply(
        a2d2inv, (np.multiply(np.multiply(a2d2inv, ymz), Ra) - np.multiply(y, Rainv))
    )
    f_015 = np.multiply(
        -a2d2inv, (np.multiply(np.multiply(a2d2inv, ymz), Ra) + np.multiply(z, Rainv))
    )
    f_115 = np.multiply(np.multiply(-a2d2inv, ypz), Rainv)
    f_215 = Rainv - np.multiply(z, f_115)
    f_125 = Rainv - np.multiply(y, f_115)

    mf = np.array([1, -1, -1, 1])
    Fintegrals = np.array(
        [
            f_003 @ mf,
            f_103 @ mf,
            f_013 @ mf,
            f_113 @ mf,
            f_213 @ mf,
            f_123 @ mf,
            f_005 @ mf,
            f_105 @ mf,
            f_015 @ mf,
            f_115 @ mf,
            f_215 @ mf,
            f_125 @ mf,
        ]
    ).T
    #                       1        2        3        4        5        6        7        8        9        10       11       12
    m4p = 0.25 * mu / np.pi
    m8p = 0.5 * m4p
    m4pn = m4p / (1 - nu)
    a2m4pn = a2 * m4pn
    a2m8p = a2 * m8p

    tdb = np.sum(np.multiply(t, b), axis=1)
    tdbv = np.array([tdb, tdb, tdb]).T
    tdbp = np.sum(np.multiply(t, bp), axis=1)
    tdbpv = np.array([tdbp, tdbp, tdbp]).T
    nddb = np.sum(np.multiply(nd, b), axis=1)
    nddbv = np.array([nddb, nddb, nddb]).T
    bxt = np.cross(b, t)
    bpxt = np.cross(bp, t)
    ndxt = np.cross(nd, t)
    bpxtdb = np.sum(np.multiply(bpxt, b), axis=1)
    bpxtdnd = np.sum(np.multiply(bpxt, nd), axis=1)
    bpxtdndv = np.array([bpxtdnd, bpxtdnd, bpxtdnd]).T
    bpxtxt = np.multiply(tdbpv, t) - bp

    I_003 = np.multiply(
        m4pn,
        (
            np.multiply(nddbv, bpxtxt)
            + np.multiply(np.array([bpxtdb, bpxtdb, bpxtdb]).T, ndxt)
            - np.multiply(bpxtdndv, bxt)
        ),
    ) - np.multiply(np.multiply(np.multiply(m4p, tdbv), tdbpv), nd)
    I_113 = np.multiply(np.multiply((m4pn - m4p), tdbv), bpxtxt)
    I_005 = (
        np.multiply(np.multiply(np.multiply(-a2m8p, tdbv), tdbpv), nd)
        - np.multiply(np.multiply(a2m4pn, bpxtdndv), bxt)
        - np.multiply(np.multiply(np.multiply(m4pn, bpxtdndv), nddbv), ndxt)
    )
    I_115 = np.multiply(np.multiply(-a2m8p, tdbv), bpxtxt) - np.multiply(
        np.multiply(np.multiply(m4pn, bpxtdndv), tdbv), ndxt
    )
    Fint_003 = Fintegrals[:, 2 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 1 - 1])
    Fint_113 = Fintegrals[:, 5 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 4 - 1])
    Fint_005 = Fintegrals[:, 8 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 7 - 1])
    Fint_115 = Fintegrals[:, 11 - 1] - np.multiply(y[:, 1 - 1], Fintegrals[:, 10 - 1])

    f4 = np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T) + np.multiply(
        I_113, np.array([Fint_113, Fint_113, Fint_113]).T
    )
    f4 = (
        f4
        + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
        + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
    )
    f4 = np.multiply(f4, np.array([oneoverL, oneoverL, oneoverL]).T)

    Fint_003 = np.multiply(y[:, 3 - 1], Fintegrals[:, 1 - 1]) - Fintegrals[:, 2 - 1]
    Fint_113 = np.multiply(y[:, 3 - 1], Fintegrals[:, 4 - 1]) - Fintegrals[:, 5 - 1]
    Fint_005 = np.multiply(y[:, 3 - 1], Fintegrals[:, 7 - 1]) - Fintegrals[:, 8 - 1]
    Fint_115 = np.multiply(y[:, 3 - 1], Fintegrals[:, 10 - 1]) - Fintegrals[:, 11 - 1]
    f3 = np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T) + np.multiply(
        I_113, np.array([Fint_113, Fint_113, Fint_113]).T
    )
    f3 = (
        f3
        + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
        + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
    )
    f3 = np.multiply(f3, np.array([oneoverL, oneoverL, oneoverL]).T)

    # Not sure if this part gets tested in the unit tests
    corindex, x1mod2, x12, x2mod2, x22, x32, x42, bp2, b2 = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(c.shape[0]):
        if np.sum(diff[i, :] ** 2) > eps * (
            np.sum(x2mod[i, :] ** 2) + np.sum(x1mod[i, :] ** 2)
        ):
            corindex.append(i)
            x1mod2.append(x1mod[i, :])
            x12.append(x1[i, :])
            x2mod2.append(x2mod[i, :])
            x22.append(x2[i, :])
            x32.append(x3[i, :])
            x42.append(x4[i, :])
            bp2.append(bp[i, :])
            b2.append(b[i, :])
    corsize = len(corindex)

    if corsize > 0:
        print("SpecialRemoteNodeForce: f3cor and f4cor")
        whocares1, whocares2, f3cor2, f4cor2 = RemoteNodeForce(
            x12, x1mod2, x32, x42, bp2, b2, a, mu, nu
        )
        whocares1, whocares2, f3cor3, f4cor3 = RemoteNodeForce(
            x2mod2, x22, x32, x42, bp2, b2, a, mu, nu
        )
        f3cor = np.zeros((c.shape[0], 3))
        f4cor = np.zeros((c.shape[0], 3))
        for i in range(corsize):
            index = corindex[i]
            f3cor[index, :] = f3cor2[i, :] + f3cor3[i, :]
            f4cor[index, :] = f4cor2[i, :] + f4cor3[i, :]
        f3 = f3 + f3cor
        f4 = f4 + f4cor

    temp = np.sum(np.multiply((x4 - x3), tp), axis=1)
    x4mod = x3 + np.multiply(np.array([temp, temp, temp]).T, tp)
    diff = x4 - x4mod
    magdiff = np.sqrt(np.sum(np.multiply(diff, diff), axis=1))
    temp = np.multiply(
        np.multiply((0.5 * cotanthetac), np.array([magdiff, magdiff, magdiff])).T, tp
    )
    x3mod = x3 + np.multiply(0.5, diff) + temp
    x4mod = x4mod + np.multiply(0.5, diff) - temp
    R = x3mod - x1
    Rdtp = np.sum(np.multiply(R, tp), axis=1)
    nd = R - np.multiply(np.array([Rdtp, Rdtp, Rdtp]).T, tp)
    d2 = np.sum(np.multiply(nd, nd), axis=1)
    r4 = np.sum(np.multiply(x4mod, tp), axis=1)
    r3 = np.sum(np.multiply(x3mod, tp), axis=1)
    s2 = np.sum(np.multiply(x2, tp), axis=1)
    s1 = np.sum(np.multiply(x1, tp), axis=1)
    y = np.array([r3, r3, r4, r4]).T
    z = np.array([-s1, -s2, -s1, -s2]).T
    a2 = a * a
    a2_d2 = a2 + d2
    temp = 1.0 / a2_d2
    a2d2inv = np.array([temp, temp, temp, temp]).T
    ypz = y + z
    ymz = y - z
    Ra = np.sqrt(np.array([a2_d2, a2_d2, a2_d2, a2_d2]).T + np.multiply(ypz, ypz))
    Rainv = 1.0 / Ra
    Log_Ra_ypz = np.log(Ra + ypz)
    f_003 = np.multiply(Ra, a2d2inv)
    f_103 = np.multiply(-0.5, (Log_Ra_ypz - np.multiply(np.multiply(ymz, Ra), a2d2inv)))
    f_013 = np.multiply(-0.5, (Log_Ra_ypz + np.multiply(np.multiply(ymz, Ra), a2d2inv)))
    f_113 = -Log_Ra_ypz
    f_213 = np.multiply(z, Log_Ra_ypz) - Ra
    f_123 = np.multiply(y, Log_Ra_ypz) - Ra
    f_005 = np.multiply(a2d2inv, (np.multiply(2.0 * a2d2inv, Ra) - Rainv))
    f_105 = np.multiply(
        a2d2inv, (np.multiply(np.multiply(a2d2inv, ymz), Ra) - np.multiply(y, Rainv))
    )
    f_015 = np.multiply(
        -a2d2inv, (np.multiply(np.multiply(a2d2inv, ymz), Ra) + np.multiply(z, Rainv))
    )
    f_115 = np.multiply(np.multiply(-a2d2inv, ypz), Rainv)
    f_215 = Rainv - np.multiply(z, f_115)
    f_125 = Rainv - np.multiply(y, f_115)
    mf = np.array([1, -1, -1, 1])
    Fintegrals = np.array(
        [
            f_003 @ mf,
            f_103 @ mf,
            f_013 @ mf,
            f_113 @ mf,
            f_213 @ mf,
            f_123 @ mf,
            f_005 @ mf,
            f_105 @ mf,
            f_015 @ mf,
            f_115 @ mf,
            f_215 @ mf,
            f_125 @ mf,
        ]
    ).T
    #                       1        2        3        4        5        6        7        8        9        10       11       12

    tpdb = np.sum(np.multiply(tp, b), axis=1)
    tpdbv = np.array([tpdb, tpdb, tpdb]).T
    tpdbp = np.sum(np.multiply(tp, bp), axis=1)
    tpdbpv = np.array([tpdbp, tpdbp, tpdbp]).T
    nddbp = np.sum(np.multiply(nd, bp), axis=1)
    nddbpv = np.array([nddbp, nddbp, nddbp]).T
    bxtp = np.cross(b, tp)
    bpxtp = np.cross(bp, tp)
    ndxtp = np.cross(nd, tp)
    bxtpdbp = np.sum(np.multiply(bxtp, bp), axis=1)
    bxtpdnd = np.sum(np.multiply(bxtp, nd), axis=1)
    bxtpdndv = np.array([bxtpdnd, bxtpdnd, bxtpdnd]).T
    bxtpxtp = np.multiply(tpdbv, tp) - b
    I_003 = np.multiply(
        m4pn,
        (
            np.multiply(nddbpv, bxtpxtp)
            + np.multiply(np.array([bxtpdbp, bxtpdbp, bxtpdbp]).T, ndxtp)
            - np.multiply(bxtpdndv, bpxtp)
        ),
    ) - np.multiply(np.multiply(np.multiply(m4p, tpdbpv), tpdbv), nd)
    I_113 = np.multiply(np.multiply((m4pn - m4p), tpdbpv), bxtpxtp)
    I_005 = (
        np.multiply(np.multiply(np.multiply(-a2m8p, tpdbpv), tpdbv), nd)
        - np.multiply(np.multiply(a2m4pn, bxtpdndv), bpxtp)
        - np.multiply(np.multiply(np.multiply(m4pn, bxtpdndv), nddbpv), ndxtp)
    )
    I_115 = np.multiply(np.multiply(-a2m8p, tpdbpv), bxtpxtp) - np.multiply(
        np.multiply(np.multiply(m4pn, bxtpdndv), tpdbpv), ndxtp
    )
    Fint_003 = Fintegrals[:, 3 - 1] - np.multiply(z[:, 1 - 1], Fintegrals[:, 1 - 1])
    Fint_113 = Fintegrals[:, 6 - 1] - np.multiply(z[:, 1 - 1], Fintegrals[:, 4 - 1])
    Fint_005 = Fintegrals[:, 9 - 1] - np.multiply(z[:, 1 - 1], Fintegrals[:, 7 - 1])
    Fint_115 = Fintegrals[:, 12 - 1] - np.multiply(z[:, 1 - 1], Fintegrals[:, 10 - 1])
    f2 = np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T) + np.multiply(
        I_113, np.array([Fint_113, Fint_113, Fint_113]).T
    )
    f2 = (
        f2
        + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
        + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
    )
    f2 = np.multiply(f2, np.array([oneoverLp, oneoverLp, oneoverLp]).T)
    Fint_003 = np.multiply(z[:, 2 - 1], Fintegrals[:, 1 - 1]) - Fintegrals[:, 3 - 1]
    Fint_113 = np.multiply(z[:, 2 - 1], Fintegrals[:, 4 - 1]) - Fintegrals[:, 6 - 1]
    Fint_005 = np.multiply(z[:, 2 - 1], Fintegrals[:, 7 - 1]) - Fintegrals[:, 9 - 1]
    Fint_115 = np.multiply(z[:, 2 - 1], Fintegrals[:, 10 - 1]) - Fintegrals[:, 12 - 1]
    f1 = np.multiply(I_003, np.array([Fint_003, Fint_003, Fint_003]).T) + np.multiply(
        I_113, np.array([Fint_113, Fint_113, Fint_113]).T
    )
    f1 = (
        f1
        + np.multiply(I_005, np.array([Fint_005, Fint_005, Fint_005]).T)
        + np.multiply(I_115, np.array([Fint_115, Fint_115, Fint_115]).T)
    )
    f1 = np.multiply(f1, np.array([oneoverLp, oneoverLp, oneoverLp]).T)

    # Not sure if this part gets tested in the unit tests
    corindex, x3mod2, x32, x4mod2, x42, x12, x22, bp2, b2 = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(c.shape[0]):
        if np.sum(diff[i, :] ** 2) > eps * (
            np.sum(x4mod[i, :] ** 2) + np.sum(x3mod[i, :] ** 2)
        ):
            corindex.append(i)
            x12.append(x1[i, :])
            x22.append(x2[i, :])
            x32.append(x3[i, :])
            x3mod2.append(x3mod[i, :])
            x42.append(x4[i, :])
            x4mod2.append(x4mod[i, :])
            bp2.append(bp[i, :])
            b2.append(b[i, :])
    corsize = len(corindex)

    if corsize > 0:
        print("SpecialRemoteNodeForce: f1cor and f2cor")
        whocares1, whocares2, f1cor2, f2cor2 = RemoteNodeForce(
            x32, x3mod2, x12, x22, b2, bp2, a, mu, nu
        )
        whocares1, whocares2, f1cor3, f2cor3 = RemoteNodeForce(
            x4mod2, x42, x12, x22, b2, bp2, a, mu, nu
        )
        f1cor = np.zeros((c.shape[0], 3))
        f2cor = np.zeros((c.shape[0], 3))
        for i in range(corsize):
            index = corindex[i]
            f1cor[index, :] = f1cor2[i, :] + f1cor3[i, :]
            f2cor[index, :] = f2cor2[i, :] + f2cor3[i, :]
        f1 = f1 + f1cor
        f2 = f2 + f2cor

    # flip back
    x1[flip, :] = x2_old[flip, :]
    x2[flip, :] = x1_old[flip, :]
    tp[flip, :] = -tp[flip, :]
    bp[flip, :] = -bp[flip, :]

    f1_old, f2_old = np.copy(f1), np.copy(f2)
    f1[flip, :] = f2_old[flip, :]
    f2[flip, :] = f1_old[flip, :]

    return f1, f2, f3, f4


################################################################################################


def python_segseg_force_vec(p1, p2, p3, p4, b1, b2, mu, nu, a):
    # change the order of input arguments: mu, nu, a -> a, mu, nu
    f1, f2, f3, f4 = RemoteNodeForce(p1, p2, p3, p4, b1, b2, a, mu, nu)
    return f1, f2, f3, f4
