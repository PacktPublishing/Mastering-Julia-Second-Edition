struct Queen
    x::Int
    y::Int
end

qhorz(qa, qb) = qa.x == qb.x;
qvert(qa, qb) = qa.y == qb.y;
qdiag(qa, qb) = abs(qa.x - qb.x) == abs(qa.y - qb.y);

qhvd(qa, qb) = qhorz(qa, qb) || qvert(qa, qb) || qdiag(qa, qb);
qany(testq, qs) = any(q -> qhvd(testq, q), qs);

function qsolve(nsqsx, nsqsy, nqs, presqs = ())
    nqs == 0 && return presqs
    for xsq in 1:nsqsx
        for ysq in 1:nsqsy
            testq = Queen(xsq, ysq)
            if !qany(testq, presqs)
                tryqs = (presqs..., testq)
                maybe = qsolve(nsqsx, nsqsy, nqs - 1, tryqs)
                maybe !== nothing && return maybe
            end
        end
    end
    return nothing
end

qsolve(nqs) = qsolve(nqs, nqs, nqs);
