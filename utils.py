import sys


class LogMessages:
    @staticmethod
    def init_log(bam: str, regions: list, outdir: str):
        return f"""Initializing UMIclusterer.
Python version: {sys.version}"
Input bam: {bam}")
Target regions: {regions}"
Output directory: {outdir}"
{"-" * 50}
"""

    @staticmethod
    def get_config(debug: bool):
        return {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "default": {
                    "format": "%(asctime)s - %(levelname)s(%(name)s): %(message)s",
                    "datefmt": "%d/%m/%Y %H:%M:%S",
                }
            },
            "handlers": {
                "file": {
                    "class": "logging.FileHandler",
                    "filename": "UMIclusterer.log",
                    "mode": "w",
                    "encoding": "utf-8",
                    "formatter": "default",
                }
            },
            "loggers": {
                "": {
                    "level": "DEBUG" if debug else "INFO",
                    "handlers": ["file"],
                }
            },
        }
