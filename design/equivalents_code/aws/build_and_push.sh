#!/bin/bash

#TAG=$(git rev-parse HEAD)
TAG=latest

REPO=equivalents

# docker build -t $REPO:$TAG .
# docker tag $REPO:$TAG 366424314563.dkr.ecr.us-east-1.amazonaws.com/$REPO:$TAG
# docker push 366424314563.dkr.ecr.us-east-1.amazonaws.com/$REPO:$TAG

CONTAINER_PROPERTIES=`cat << EOF
{
	"image": "366424314563.dkr.ecr.us-east-1.amazonaws.com/$REPO:$TAG",
	"vcpus": 1,
	"memory": 1024,
	"environment": [
		{
			"name": "D_MAX",
			"value": "100"
		}
	],
	"jobRoleArn": "arn:aws:iam::366424314563:role/double-descent-batch-s3"
}
EOF`
CMD="aws batch register-job-definition \
	--job-definition-name $REPO \
	--type container \
	--container-properties '$CONTAINER_PROPERTIES'"
echo $CMD
RESULT=$(eval $CMD)
JOB_ARN=$(echo $RESULT | jq '.jobDefinitionArn' -r)

SIZE=$(python run_experiment.py)
aws batch submit-job \
	--job-name $REPO \
	--job-queue default \
	--array-properties size=$SIZE \
	--job-definition $JOB_ARN
