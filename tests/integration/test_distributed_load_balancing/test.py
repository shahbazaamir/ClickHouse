# pylint: disable=unused-argument
# pylint: disable=redefined-outer-name
# pylint: disable=line-too-long

import uuid

import pytest

from helpers.cluster import ClickHouseCluster

cluster = ClickHouseCluster(__file__)

n1 = cluster.add_instance(
    "n1",
    main_configs=["configs/remote_servers.xml"],
    user_configs=["configs/users.xml"],
)
n2 = cluster.add_instance(
    "n2",
    main_configs=["configs/remote_servers.xml"],
    user_configs=["configs/users.xml"],
)
n3 = cluster.add_instance(
    "n3",
    main_configs=["configs/remote_servers.xml"],
    user_configs=["configs/users.xml"],
)

nodes = len(cluster.instances)
queries = nodes * 10


# SYSTEM RELOAD CONFIG will reset some attributes of the nodes in cluster
# - error_count
# - last_used (round_robing)
#
# This is required to avoid interference results of one test to another
@pytest.fixture(scope="function", autouse=True)
def test_setup():
    for n in list(cluster.instances.values()):
        n.query("SYSTEM RELOAD CONFIG")


def bootstrap():
    for n in list(cluster.instances.values()):
        n.query("DROP TABLE IF EXISTS data")
        n.query("DROP TABLE IF EXISTS dist")
        n.query("CREATE TABLE data (key Int) Engine=Memory()")
        n.query(
            """
        CREATE TABLE dist AS data
        Engine=Distributed(
            replicas_cluster,
            currentDatabase(),
            data)
        """
        )
        n.query(
            """
        CREATE TABLE dist_priority AS data
        Engine=Distributed(
            replicas_priority_cluster,
            currentDatabase(),
            data)
        """
        )
        n.query(
            """
        CREATE TABLE dist_priority_negative AS data
        Engine=Distributed(
            replicas_priority_negative_cluster,
            currentDatabase(),
            data)
        """
        )


def make_uuid():
    return uuid.uuid4().hex


@pytest.fixture(scope="module", autouse=True)
def start_cluster():
    try:
        cluster.start()
        bootstrap()
        yield cluster
    finally:
        cluster.shutdown()


def get_node(query_node, table="dist", *args, **kwargs):
    query_id = make_uuid()

    settings = {
        "query_id": query_id,
        "log_queries": 1,
        "log_queries_min_type": "QUERY_START",
        "prefer_localhost_replica": 0,
        "max_parallel_replicas": 1,
    }
    if "settings" not in kwargs:
        kwargs["settings"] = settings
    else:
        kwargs["settings"].update(settings)

    query_node.query("SELECT * FROM " + table, *args, **kwargs)

    for n in list(cluster.instances.values()):
        n.query("SYSTEM FLUSH LOGS")

    rows = query_node.query(
        """
    SELECT hostName()
    FROM cluster(shards_cluster, system.query_log)
    WHERE
        initial_query_id = '{query_id}' AND
        is_initial_query = 0 AND
        type = 'QueryFinish'
    ORDER BY event_date DESC, event_time DESC
    LIMIT 1
    """.format(
            query_id=query_id
        )
    )
    return rows.strip()


# TODO: right now random distribution looks bad, but works
def test_load_balancing_default():
    unique_nodes = set()
    for _ in range(0, queries):
        unique_nodes.add(get_node(n1, settings={"load_balancing": "random"}))
    assert len(unique_nodes) == nodes, unique_nodes


def test_load_balancing_nearest_hostname():
    unique_nodes = set()
    for _ in range(0, queries):
        unique_nodes.add(get_node(n1, settings={"load_balancing": "nearest_hostname"}))
    assert len(unique_nodes) == 1, unique_nodes
    assert unique_nodes == set(["n1"])


def test_load_balancing_hostname_levenshtein_distance():
    unique_nodes = set()
    for _ in range(0, queries):
        unique_nodes.add(
            get_node(n1, settings={"load_balancing": "hostname_levenshtein_distance"})
        )
    assert len(unique_nodes) == 1, unique_nodes
    assert unique_nodes == set(["n1"])


def test_load_balancing_in_order():
    unique_nodes = set()
    for _ in range(0, queries):
        unique_nodes.add(get_node(n1, settings={"load_balancing": "in_order"}))
    assert len(unique_nodes) == 1, unique_nodes
    assert unique_nodes == set(["n1"])


def test_load_balancing_first_or_random():
    unique_nodes = set()
    for _ in range(0, queries):
        unique_nodes.add(get_node(n1, settings={"load_balancing": "first_or_random"}))
    assert len(unique_nodes) == 1, unique_nodes
    assert unique_nodes == set(["n1"])


def test_load_balancing_round_robin():
    unique_nodes = set()
    for _ in range(0, nodes):
        unique_nodes.add(get_node(n1, settings={"load_balancing": "round_robin"}))
    assert len(unique_nodes) == nodes, unique_nodes
    assert unique_nodes == set(["n1", "n2", "n3"])


@pytest.mark.parametrize(
    "dist_table",
    [
        ("dist_priority"),
        ("dist_priority_negative"),
    ],
)
def test_load_balancing_priority_round_robin(dist_table):
    unique_nodes = set()
    for _ in range(0, nodes):
        unique_nodes.add(
            get_node(n1, dist_table, settings={"load_balancing": "round_robin"})
        )
    assert len(unique_nodes) == 2, unique_nodes
    # n2 has bigger priority in config
    assert unique_nodes == set(["n1", "n3"])


def test_distributed_replica_max_ignored_errors():
    settings = {
        "use_hedged_requests": 0,
        "load_balancing": "in_order",
        "prefer_localhost_replica": 0,
        "connect_timeout": 2,
        "receive_timeout": 2,
        "send_timeout": 2,
        "tcp_keep_alive_timeout": 2,
        "distributed_replica_max_ignored_errors": 0,
        "distributed_replica_error_half_life": 60,
        "max_parallel_replicas": 1,
    }

    # initiate connection (if started only this test)
    n2.query("SELECT * FROM dist", settings=settings)

    with cluster.pause_container("n1"):
        # n1 paused -- skipping, and increment error_count for n1
        # but the query succeeds, no need in query_and_get_error()
        n2.query("SELECT * FROM dist", settings=settings)
        # XXX: due to config reloading we need second time (sigh)
        n2.query("SELECT * FROM dist", settings=settings)
        # check error_count for n1
        assert (
            int(
                n2.query(
                    """
        SELECT errors_count FROM system.clusters
        WHERE cluster = 'replicas_cluster' AND host_name = 'n1'
        """,
                    settings=settings,
                )
            )
            == 1
        )

    # still n2
    assert get_node(n2, settings=settings) == "n2"
    # now n1
    settings["distributed_replica_max_ignored_errors"] = 1
    assert get_node(n2, settings=settings) == "n1"
