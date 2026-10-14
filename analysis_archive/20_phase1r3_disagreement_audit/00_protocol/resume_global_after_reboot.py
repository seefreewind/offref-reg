"""Resume only the unfinished global worker, then run guarded completion steps."""
from pathlib import Path
import datetime
import json
import subprocess
import sys

R = Path(__file__).resolve().parents[1]
PY = sys.executable
sys.path.insert(0, str(R / '00_protocol'))
from retry_transient_failures import archive_for_one_retry


def state(status, **details):
    (R / '14_reports/EXECUTION_STATUS.json').write_text(json.dumps({
        'status': status, 'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'Phase2': 'NOT_ENTERED', **details,
    }, indent=2))


def main():
    folders = {'GLOBAL_CORE': '02_track_harmonization', 'ENSEMBLE': '03_within_model_ensembles',
               'PLACEMENT': '04_context_robustness', 'CONTEXT': '04_context_robustness',
               'GNRS_GENE_CENTER': '08_truth_coverage', 'REFERENCE': '07_reference_benchmark'}
    jobs = json.loads((R / '00_protocol/PREDICTION_JOBS_FROZEN.json').read_text())
    missing_other = [j['audit_job_ID'] for j in jobs
                     if not (j['model'] == 'AlphaGenome' and j['stage'] == 'GLOBAL_CORE')
                     and not (R / folders[j['stage']] / j['model'] / (j['audit_job_ID'] + '.json')).exists()]
    assert not missing_other, 'NON_GLOBAL_JOBS_MISSING'
    for attempt in range(2):
        if attempt:
            archive_for_one_retry('AlphaGenome', ['GLOBAL_CORE'], 1, 0)
        with (R / '01_model_metadata/AG_GLOBAL_post_reboot.log').open('a') as log:
            worker = subprocess.Popen([PY, str(R / '00_protocol/run_predictions.py'),
                                       '--model', 'AlphaGenome', '--stages', 'GLOBAL_CORE'],
                                      stdout=log, stderr=subprocess.STDOUT)
            state('FORMAL_PREDICTIONS_RUNNING', supervisor_PID=__import__('os').getpid(),
                  global_PID=worker.pid, resumption_attempt=attempt,
                  reason='Host reboot; completed non-global batches preserved')
            worker.wait()
            assert worker.returncode == 0, 'GLOBAL_WORKER_PROCESS_FAILED'
    state('PREDICTIONS_COMPLETE_ANALYSIS_PENDING')
    for script in ['11_statistics/analyze.py', '09_disagreement_models/fit_models.py',
                   '11_statistics/finalize_endpoints.py', '12_figures/render_figures.py',
                   '14_reports/build_report.py']:
        with (R / '14_reports/post_reboot_completion.log').open('a') as log:
            subprocess.run([PY, str(R / script)], stdout=log, stderr=subprocess.STDOUT, check=True)
    state('FIGURES_AND_DRAFT_REPORT_COMPLETE_MANUAL_QA_REQUIRED')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        state('CONTROLLER_FAILED', error_type=type(error).__name__)
        raise
