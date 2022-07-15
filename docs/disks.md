# Additional disk support

When we create our EC2 instances we configured two additional drives. These get initialized and mounted using the following commmands. The `-E discard` enables SSD TRIM operations, `-m 1` allocates 1% of the available space to the super-user, and `-O ^has_journal` disables journaling. We disable journaling because these are effectively short-lived EC2 instances and if they crash we won't try to recover them.

Note that the following does not add the drives to `/etc/fstab` which means that it's your responsibility to re-mount them in the case your EC2 restarts.

```bash
sudo mkfs.ext4 -E discard -m 1 -O ^has_journal /dev/nvme1n1
sudo mkfs.ext4 -E discard -m 1 -O ^has_journal /dev/nvme2n1

sudo mount /dev/nvme1n1 ${HOME}/pipeline
sudo mount /dev/nvme2n1 ${HOME}/stats

sudo chown -R ubuntu:ubuntu ${HOME}/pipeline
sudo chown -R ubuntu:ubuntu ${HOME}/stats
```

Setting up additional swap space happens on initial EC2 boot. You will need to run `swapon` again if the instance is rebooted or a snapshot is taken (which disconnects drives).

```bash
sudo dd if=/dev/zero of=${HOME}/stats/swap bs=1M count=131072
sudo chmod 0600 ${HOME}/stats/swap
sudo mkswap ${HOME}/stats/swap
sudo swapon ${HOME}/stats/swap
```