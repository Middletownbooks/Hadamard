r"""The exact joint test of one y_B under a hard memory cap [M].

nd_joint_exhaust.JointTest.test can need more memory than the machine has for some y_B (its layer-option and pair
enumerations grow with the y_B; y_B #11 of 495 x83 was killed at 8.6 GB with the pair lists capped at 4e6).  Here the
test runs in a forked child whose address space is capped (RLIMIT_AS); the child returns (rec, sols) through a pipe.
A MemoryError, any other exception, or a dead child is reported, never an out-of-memory kill of the calling run.
capped_test(X, yb, maxsol, maxpairs, mem_gb) -> (status, rec, sols), status in {'ok', 'memory', 'error', 'died'}."""
import multiprocessing as mp, resource, traceback

def _child(conn, X, yb, maxsol, maxpairs, mem_bytes):
    try:
        resource.setrlimit(resource.RLIMIT_AS, (mem_bytes, mem_bytes))
        rec, sols = X.J.test(yb, maxsol=maxsol, maxpairs=maxpairs)
        conn.send(('ok', rec, sols))
    except MemoryError:
        conn.send(('memory', None, None))
    except BaseException:
        try: conn.send(('error', traceback.format_exc()[-400:], None))
        except BaseException: pass
    finally:
        conn.close()

def capped_test(X, yb, maxsol=10 ** 7, maxpairs=4 * 10 ** 6, mem_gb=6.0):
    ctx = mp.get_context('fork')
    rd, wr = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_child, args=(wr, X, yb, maxsol, maxpairs, int(mem_gb * 2 ** 30)))
    p.start(); wr.close()
    try:
        status, rec, sols = rd.recv()
    except (EOFError, OSError):
        status, rec, sols = 'died', None, None
    p.join()
    if status == 'ok': return status, rec, sols
    return status, dict(joint=None, note=f"joint test {status} (cap {mem_gb} GB address space)" +
                        (f": {rec}" if status == 'error' else "")), []
