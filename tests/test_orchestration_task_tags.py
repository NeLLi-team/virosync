from virosync.orchestration._flows.single_genome import run_single_genome_task
from virosync.orchestration.tasks import (
    classify_jelly_roll_task,
    generate_outputs_task,
    generate_proteome_task,
    gene_taxonomy_batch_task,
    hhg_seeding_task,
    interproscan_batch_task,
    marker_validation_task,
    mask_genome_task,
    region_assembly_task,
    taxonomy_expansion_task,
    verify_eve_candidates_batched_task,
    verify_eve_task,
)


def test_orchestration_tasks_are_plain_python_callables() -> None:
    tasks = [
        mask_genome_task,
        generate_proteome_task,
        hhg_seeding_task,
        marker_validation_task,
        region_assembly_task,
        taxonomy_expansion_task,
        gene_taxonomy_batch_task,
        interproscan_batch_task,
        verify_eve_task,
        verify_eve_candidates_batched_task,
        classify_jelly_roll_task,
        generate_outputs_task,
        run_single_genome_task,
    ]

    for task_func in tasks:
        assert callable(task_func)
        assert not hasattr(task_func, "submit")
        assert not hasattr(task_func, "tags")
