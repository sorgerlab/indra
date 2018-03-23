import boto3
from indra.db.util import get_db_statistics
from datetime import datetime

def main():
    utcnow = datetime.utcnow()
    fname = "Primary_Database_Status_Report_%s.txt" % utcnow.strftime("%Y%m%d")
    print("Creating report in: %s." % fname)

    print("\nBegin Report============\n")
    get_db_statistics(fname)
    print("\nEnd Report==============\n")

    print("Saving record to s3.")
    s3 = boto3.client('s3')
    with open(fname, 'rb') as f:
        s3.put_object(Body=f, Bucket='bigmech', Key='indra-db/reports/%s' % fname)
    return

if __name__ == '__main__':
    main()
