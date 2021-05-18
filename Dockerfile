FROM 579831337053.dkr.ecr.us-west-2.amazonaws.com/irx-r-base:6

# Show all lines of testing output
ENV _R_CHECK_TESTS_NLINES_=0

COPY ./ /src
WORKDIR /src

# Install dependencies
RUN R CMD INSTALL ./PKPDsim2 --no-multiarch --no-docs
RUN rm -rf ./PKPDsim2
