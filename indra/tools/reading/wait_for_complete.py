import sys
from indra.tools.reading.submit_reading_pipeline_aws import wait_for_complete

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Missing argument. Must specify queue.")
        print("Usage: python wait_for_complete.py QUEUE")
        sys.exit(1)
    wait_for_complete(sys.argv[-1])
